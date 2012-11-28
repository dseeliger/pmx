#include "Python.h"
#include "arrayobject.h"
#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h" //Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
#include "NGMX_XTC.h"




/* #### numeric Python interface  #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef NGMX_XTCMethods[] = {
    {"xdrfile_NT", 		xdrfile_NT, 	METH_VARARGS},
    {"xdrfile_read", 	xdrfile_read, 	METH_VARARGS},
    //{"xdrfile_readI", 	xdrfile_readI, 	METH_VARARGS},
    {"xdrfile_write", 	xdrfile_write, 	METH_VARARGS},
    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

/* ==== Initialize the C_test functions ====================== */
// Module name must be _C_arraytest in compile and linked
void initNGMX_XTC()  {
    (void) Py_InitModule("NGMX_XTC", NGMX_XTCMethods);
    import_array();  // Must be present for NumPy.  Called first after above line.
}


struct XDRFILE
{
    FILE *   fp;       //< pointer to standard C library file handle /
    void *   xdr;      //< pointer to corresponding XDR handle --> not used here!  /
    char     mode;     //< r=read, w=write, a=append                 /
    int *    buf1;     //< Buffer for internal use                   /
    int      buf1size; //< Current allocated length of buf1          /
    int *    buf2;     //< Buffer for internal use                   /
    int      buf2size; //< Current allocated length of buf2          /
};


struct XDR_HEAD
{
	//23*4 = 92 bytes long
	int 			Magic;		//Magic number should be 1995 (always)
	int 			N1;			//Number of Atoms
	int 			frame;		//Frame number
	float			time;		//time for the frame
	float			box[9];		//tensor for the box dimensions
	int				N2;			//Number of Atoms repeated
	float			precision;	//the XTC precision (as stated in the mdp file)
	int				minint[3];	//minimum integer for each dimension
	int				maxint[3];	//maximum integer for each dimension
	int				smallidx;	// ????
	unsigned int	data_len;	//length (in byte) of the (compressed) databuffer "j += (j%4) ? 0 : 4-(j%4)" --> to avoid frame shifts if no multiple of 4 (32bit)
};



/* #### Vector Extensions ############################## */

/* ==== vector function - manipulate vector in place ======================*/


//Stolen from xdrfile.c ////////////////////

static const int magicints[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031,
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216
};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

static int
sizeofints(int num_of_ints, unsigned int sizes[])
{
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++)
    {
		tmp = 0;
		for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0)
        {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num)
    {
		num_of_bits++;
		num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}

static int
sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;

    while (size >= num && num_of_bits < 32)
    {
		num_of_bits++;
		num <<= 1;
    }
    return num_of_bits;
}

static void
encodebits(int buf[], int num_of_bits, int num)
{

    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8)
    {
		lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
		cbuf[cnt++] = lastbyte >> lastbits;
		num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
		lastbyte = (lastbyte << num_of_bits) | num;
		lastbits += num_of_bits;
		if (lastbits >= 8)
        {
			lastbits -= 8;
			cbuf[cnt++] = lastbyte >> lastbits;
		}
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0)
    {
		cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

static void
encodeints(int buf[], int num_of_ints, int num_of_bits,
		   unsigned int sizes[], unsigned int nums[])
{

    int i;
    unsigned int bytes[32], num_of_bytes, bytecnt, tmp;

    tmp = nums[0];
    num_of_bytes = 0;
    do
    {
		bytes[num_of_bytes++] = tmp & 0xff;
		tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++)
    {
		if (nums[i] >= sizes[i])
        {
			fprintf(stderr,"major breakdown in encodeints - num %u doesn't "
					"match size %u\n", nums[i], sizes[i]);
			abort();
		}
		/* use one step multiply */
		tmp = nums[i];
		for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0)
        {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8)
    {
		for (i = 0; i < num_of_bytes; i++)
        {
			encodebits(buf, 8, bytes[i]);
		}
		encodebits(buf, num_of_bits - num_of_bytes * 8, 0);
    }
    else
    {
		for (i = 0; i < num_of_bytes-1; i++)
        {
			encodebits(buf, 8, bytes[i]);
		}
		encodebits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}

static int
decodebits(int buf[], int num_of_bits)
{

    int cnt, num;
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];

    num = 0;
    while (num_of_bits >= 8)
    {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
		num_of_bits -=8;
    }
    if (num_of_bits > 0)
    {
		if (lastbits < num_of_bits)
        {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
		}
		lastbits -= num_of_bits;
		num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}


static void
decodeints(int buf[], int num_of_ints, int num_of_bits,
		   unsigned int sizes[], int nums[])
{

	int bytes[32];
	int i, j, num_of_bytes, p, num;

	bytes[1] = bytes[2] = bytes[3] = 0;
	num_of_bytes = 0;
	while (num_of_bits > 8)
    {
		bytes[num_of_bytes++] = decodebits(buf, 8);
		num_of_bits -= 8;
	}
	if (num_of_bits > 0)
    {
		bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
	}
	for (i = num_of_ints-1; i > 0; i--)
    {
		num = 0;
		for (j = num_of_bytes-1; j >=0; j--)
        {
			num = (num << 8) | bytes[j];
			p = num / sizes[i];
			bytes[j] = p;
			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}


////////////end of stolen Part


/////Function for the Python routines

int skip_frame(XDRFILE* xfp, XDR_HEAD* header, long file_byte_len)
{
	int dl;
	if (header->N1 < 9) //if less then 9 atoms
	{
		dl = header->N1 * 4;
		if (dl + ftell(xfp->fp) > file_byte_len)
		{
			PyErr_SetString(PyExc_EOFError, "File not long enough to read another set of data!");
			return NULL;
		}
			return NULL;
		fseek(xfp->fp, dl, SEEK_CUR);
		return 1;
	}
	//more than 9 atoms
	dl = header->data_len;
	if( dl % 4 ) dl += 4 - (dl % 4);
	if (dl + ftell(xfp->fp) > file_byte_len)
		return NULL;
	fseek(xfp->fp, dl, SEEK_CUR);
	return 1;

}

int read_header(XDRFILE* xfp, XDR_HEAD* header, long file_byte_len)
{
	if (  (ftell(xfp->fp) + 23*4) > file_byte_len  )
	{
		PyErr_SetString(PyExc_EOFError, "File not long enough to read another header!");
		return NULL;
	};

	int	*i_header;
	i_header = header;	//"cast" to pointer of integer to make it compatible with the reading routine
	if ( !xdrfile_read_int(i_header, 23, xfp) )
	{
		PyErr_SetString(PyExc_MemoryError, "Unable to read header");
		return NULL;
	};


	if( header->Magic != 1995)
	{
		PyErr_SetString(PyExc_ValueError, "Magic number Error!");
		return NULL; //check magic number
	};

	return 1; //All clear!
}

static PyObject *xdrfile_NT(PyObject *self, PyObject *args)
{
	const char *file_path;
	XDRFILE* xfp;
	XDR_HEAD *header;
	int T=0, N=-1;
	long file_byte_len;
	if (!PyArg_ParseTuple(args, "s", &file_path))  return 0;

	xfp = xdrfile_open(file_path, "r");
	// storing the length of the file in bytes to determine when we are at the end
	fseek(xfp->fp, 0L, SEEK_END);
	file_byte_len = ftell(xfp->fp);
	fseek(xfp->fp, 0L, SEEK_SET);

	//making space for the header
	header = (XDR_HEAD*) malloc( sizeof(XDR_HEAD) );

	// going through the File
	while( 1 )
	{

		//read header
		if (  !read_header( xfp, header, file_byte_len )  )
		{
			if ( PyErr_ExceptionMatches(PyExc_EOFError) ) //Catch error EOF
			{
				PyErr_Clear();
				break;
			}
			free(header);		//If error occured leave with the error
			xdrfile_close(xfp); //but clean up!
			return NULL;
		}

		//see if number stays constant else ERROR
		if(N == -1)		//only first frame...
			N = header->N1;
		if(header->N1    != N )
		{
			PyErr_SetString(PyExc_ValueError, "Different numbers of atoms in two different frames");
			free(header);
			xdrfile_close(xfp);
			return NULL;
		}

		//skip over data
		if (  !skip_frame( xfp, header, file_byte_len )  )
		{
			if ( PyErr_ExceptionMatches(PyExc_EOFError) ) //Catch error EOF
			{
				PyErr_Clear();
				break;
			}
			free(header);		//If error occured leave with the error
			xdrfile_close(xfp); //but clean up!
			return NULL;
		}
		T++;
	}
	free(header);
	xdrfile_close(xfp);
	return Py_BuildValue("i,i", N, T);
}

static PyObject *xdrfile_read(PyObject *self, PyObject *args)
{
//self.filename, DatObj.x, DatObj.time, DatObj.BOX, start, stop, skip, verbose=True

	const char *file_path;				//file path of XTC file
	long file_byte_len, temp_file_pos;	//the length of the file in bytes
	XDRFILE* xfp;						//file handle
	XDR_HEAD* header;					//file header (with all important info)

	PyArrayObject *x, *time, *box, *ndx, *step;//Python arrays
	float *time_data, *box_data;    //pointer to the data of the python arrays
	unsigned long long *step_data, *ndx_ptr;
	int start, stop, skip, verb=0;		//where to start end and how many frames to skip (and if verbose or not)
	unsigned int N_stride, ndx_i, ndx_len;		//length
	float scaling=1.0;					//scales all coordinates by this number (useful if you want to go from nm->A)

	int T=0, Tr=0;						//Frame counter looping variable; Read frame to address the arrays

	//// Standard lib initializations
	int *minint, *maxint;
	int smallidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int k, *buf2, lsize, flag;
	int smallnum, smaller, i, is_smaller, run;
	float *lfp, inv_precision;
	int tmp, thiscoord[3],  prevcoord[3];
	unsigned int bitsize;
	////////////////

	if (!PyArg_ParseTuple(args, "sO!O!O!O!O!iii|fi:xdrfile_read",
    	&file_path,
        &PyArray_Type, &x,
        &PyArray_Type, &box,
        &PyArray_Type, &time,
        &PyArray_Type, &step,
        &PyArray_Type, &ndx,
        &start, &stop, &skip,
        &scaling, &verb				))  return Py_BuildValue("i", 0);

	//opening file handle...
	xfp = xdrfile_open(file_path, "r");
	//...and getting length of the file
	fseek(xfp->fp, 0L, SEEK_END);
	file_byte_len = ftell(xfp->fp);
	fseek(xfp->fp, 0L, SEEK_SET);

	//making space for the header
	header = (XDR_HEAD*) malloc( sizeof(XDR_HEAD) );

	//getting access to the Data of the arrays
	time_data = time->data;
	box_data  = box->data;
	step_data = step->data;
	ndx_len   = PyArray_DIM(ndx, 0); //length of the index vector
	N_stride 	= PyArray_STRIDE(x, 1) / 4; // not sure why we divide by 4, but that way it works... //&&&



	///Initializations independent of time frame or Atom number (we assume that the number of coordinates won't change with the time frame!)
	if( !read_header(xfp, header, file_byte_len) )
		return NULL;				//read the header once to initialize variables
	fseek(xfp->fp, 0L, SEEK_SET);	//go back to the start..

	lsize = header->N1;							//Number of Atoms
	size3 = lsize * 3;							//Number of Floats needed (#Atoms * 3dim)
	inv_precision = scaling / header->precision;


	//Buffer2 for encoded data
	if(  (xfp->buf2=(int *)malloc(sizeof(int)*size3*1.2)) == NULL  )
	        {
				PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory for decompressing coordinates.\n");
				xdrfile_close(xfp);
				free(header);
				return NULL;
			}
	buf2=xfp->buf2;


	while(T < start) //skip until we are at the frame we want to start reading from to the start position
	{
		//skip over data
		if( !read_header(xfp, header, file_byte_len) )
			return NULL;
		if( !skip_frame(xfp, header, file_byte_len) )
			return NULL;
		T++;
	}


	while( T < stop ) //now read every go on with the looping but read in every "skip"th frame
	{

		//		Read HEADER
		if( !read_header(xfp, header, file_byte_len) )
		{
			xdrfile_close(xfp); //clean exit ...
			free(header);
			//PyErr_SetString(PyExc_MemoryError, "could not read header");
			return NULL;
		}

		//		TEST FOR SKIP AND DO SO
		if( (T-start) % skip) //if expr==0 --> False
		{
			//skip over data
			if( !skip_frame(xfp, header, file_byte_len) )
			{
				PyErr_SetString(PyExc_MemoryError, "could not jump");
				xdrfile_close(xfp);
				free(header);
				return NULL;
			}
			T++;
			continue; //avoid the else clause...
		}


		// IF WE ARE HERE IT'S READING TIME
		if(verb)
			printf("reading frame %i\r",T);

		// write box, and time info to the array
		time_data[Tr] = header->time;
		step_data[Tr] = header->frame;
		for(i=0; i<9; i++)
			box_data[ Tr*9 + i ] = header->box[i] * scaling;


	    //TODO Check if number of atoms stays constant ...
		//TODO Decompression for less than 9 Atoms



		temp_file_pos	= ftell(xfp->fp);
		bitsizeint[0] 	= bitsizeint[1] =  bitsizeint[2] = 0;
		buf2[0] 		= buf2[1] = buf2[2] = 0;
		maxint 			= header->maxint;
		minint 			= header->minint;
		sizeint[0] 		= maxint[0] - minint[0]+1;
		sizeint[1] 		= maxint[1] - minint[1]+1;
		sizeint[2] 		= maxint[2] - minint[2]+1;
		lfp 			= PyArray_GETPTR1(x, Tr);		//pointer to the right Time frame of coord array //&&&
		ndx_ptr  		= ndx->data;						//pointer to the current index that we want to put to the output
		ndx_i	  		= 0;	//counter for the next index we want to actually write out
		run 			= 0;	//counter for runs (multiple decodes)
		i 				= 0;	//counter for "runs" (multiple coordinates) and the atom counter


		/* check if one of the sizes is to big to be multiplied */
		if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) //never seen it happen...
	    {
			bitsizeint[0] = sizeofint(sizeint[0]);
			bitsizeint[1] = sizeofint(sizeint[1]);
			bitsizeint[2] = sizeofint(sizeint[2]);
			bitsize = 0; /* flag the use of large sizes */
		}

	    else
			bitsize = sizeofints(3, sizeint);


		smallidx = header->smallidx;
		tmp = (FIRSTIDX>smallidx-1) ? FIRSTIDX : smallidx-1;
		smaller = magicints[tmp] / 2;
		smallnum = magicints[smallidx] / 2;
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;


		//read the compressed data from the file into buffer 2
		if (  xdrfile_read_opaque((char *)&(buf2[3]),  (unsigned int)header->data_len,  xfp) == 0  )
		{
			PyErr_SetString(PyExc_BufferError, "Unable to read compressed data.\n");
			xdrfile_close(xfp);
			free(header);
			return NULL;
		}

		while(ndx_i < ndx_len) //only go until
	    {

			if (bitsize == 0)
	        {
				thiscoord[0] = decodebits(buf2, bitsizeint[0]);
				thiscoord[1] = decodebits(buf2, bitsizeint[1]);
				thiscoord[2] = decodebits(buf2, bitsizeint[2]);
			}
	        else
	        {
				decodeints(buf2, 3, bitsize, sizeint, thiscoord);
			}

			thiscoord[0] += minint[0];
			thiscoord[1] += minint[1];
			thiscoord[2] += minint[2];

			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];


			flag = decodebits(buf2, 1);
			is_smaller = 0;
			if (flag == 1)
	        {
				run = decodebits(buf2, 5);
				is_smaller = run % 3;
				run -= is_smaller;
				is_smaller--;
			}


			if (run > 0)
	        {
				//thiscoord += 3;
				for (k = 0; k < run; k+=3)
	            {
					decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
					thiscoord[0] += prevcoord[0] - smallnum;
					thiscoord[1] += prevcoord[1] - smallnum;
					thiscoord[2] += prevcoord[2] - smallnum;
					if (k == 0) {
						// interchange first with second atom for better
						// compression of water molecules

						tmp = thiscoord[0]; thiscoord[0] = prevcoord[0]; prevcoord[0] = tmp;
						tmp = thiscoord[1]; thiscoord[1] = prevcoord[1]; prevcoord[1] = tmp;
						tmp = thiscoord[2]; thiscoord[2] = prevcoord[2]; prevcoord[2] = tmp;


						if( ndx_ptr[ndx_i] == i ) //only write to the array if we want it!
						{
							lfp[0] = prevcoord[0] * inv_precision;	//write to the array
							lfp[1] = prevcoord[1] * inv_precision;
							lfp[2] = prevcoord[2] * inv_precision;
							lfp += N_stride; 			//go to the next vector position in the array
							ndx_i++;								//go to the next index we want to read
						}
						i++;										//always go to the next Atom

					}
					else
					{
						prevcoord[0] = thiscoord[0];
						prevcoord[1] = thiscoord[1];
						prevcoord[2] = thiscoord[2];
					}
					if( ndx_ptr[ndx_i] == i )  //only write to the array if we want it!
					{
						lfp[0] = thiscoord[0] * inv_precision;
						lfp[1] = thiscoord[1] * inv_precision;
						lfp[2] = thiscoord[2] * inv_precision;
						lfp += N_stride;
						ndx_i++;
					}
					i++;//end writing
				}
			}
	        else
	        {

				if( ndx_ptr[ndx_i] == i )  //only write to the array if we want it!
				{
					lfp[0] = thiscoord[0] * inv_precision;
					lfp[1] = thiscoord[1] * inv_precision;
					lfp[2] = thiscoord[2] * inv_precision;
					lfp += N_stride;
					ndx_i++;
				}
				i++;//end writing

			}
			smallidx += is_smaller;
			if (is_smaller < 0)
	        {
				smallnum = smaller;

				if (smallidx > FIRSTIDX)
	            {
					smaller = magicints[smallidx - 1] /2;
				}
	            else
	            {
					smaller = 0;
				}
			}
	        else if (is_smaller > 0)
	        {
				smaller = smallnum;
				smallnum = magicints[smallidx] / 2;
			}
			sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
		}
		fseek(xfp->fp, temp_file_pos , SEEK_SET); //jump to the start of this frame (in case we did not read the entire frame and skipped over the last atoms)
		skip_frame(xfp, header, file_byte_len);	  // and skip over it...
		T++; Tr++; //increment the "frame number" and the "read frames number"
	}
	printf("\nClosing....\n");
	xdrfile_close(xfp);
	free(header);
	return Py_BuildValue("i", 1);
};


static PyObject *xdrfile_write(PyObject *self, PyObject *args)
{
	const char 			*file_path;						//file path of XTC file
	XDRFILE 			*xfp;							//file handle
	XDR_HEAD 			*header;						//file header (with all important info)
	int   				head_int[3];					//variables to write the headder (magic, natoms, step)
	float 				head_float[10];					//variables to write the headder (time, 3x3 box )
	PyArrayObject		*x, *time, *box, *step;				//Python arrays
	float 				*x_data, *time_data, *box_data;	//pointers to tha data of the Python array
	int					X_Tstride, X_Nstride, X_Cstride; //the strides to move through the array
	long int			*step_data;
	float 				scaling=10.0, start_prec=10000.0;
	int					start, stop, skip, verb;
	unsigned int	   	T_frames, N_atoms;
	long unsigned int  	T=0,  N=0;


	if (!PyArg_ParseTuple(args, "sO!O!O!O!|fif:xdrfile_write",
    	&file_path,
        &PyArray_Type, &x,
        &PyArray_Type, &box,
        &PyArray_Type, &time,
        &PyArray_Type, &step,
        &scaling, &verb, &start_prec)			)  return Py_BuildValue("i", 0);

	//get the right handles on the python arrays
	T_frames 	= (int) PyArray_DIM(x, 0); //number of time frames
	N_atoms		= (int) PyArray_DIM(x, 1); //number of atoms
	time_data	= PyArray_DATA( time );
	box_data    = PyArray_DATA( box );
	step_data   = PyArray_DATA( step );
	X_Tstride   = PyArray_STRIDE(x, 0)/sizeof(float); //stride along time frame axis
	X_Nstride   = PyArray_STRIDE(x, 1)/sizeof(float); //       along atoms axis
	X_Cstride   = PyArray_STRIDE(x, 2)/sizeof(float); //       along coordinates axis (X,Y,Z)


	xfp = xdrfile_open(file_path, "w+"); //open file to write...

	//definitions from xdrfile.c
	int minint[3], maxint[3], mindiff, *lip, diff;
	int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
	int minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
	int k, *buf1, *buf2;
	int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
	float *ptr, lf, precision;
	int tmp, tmpsum, *thiscoord,  prevcoord[3];
	unsigned int tmpcoord[30];
	int errval=1;
	unsigned int bitsize;

	//hope to get rid of
	int size;
	float *lfp;

	if(verb)
		printf("Writing XTC file \"%s\":\n", file_path);


	// From the xdrfile.c routine that can be done outside the T loop
	size3 = N_atoms*3;

	head_int[0] = 1995;
	head_int[1] = N_atoms;

	if(size3>xfp->buf1size)
    {
		if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL)
        {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL)
        {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
	}

	for(T=0; T<T_frames; T++) //loop through timeframes
	{

		//Write the header
		head_int[2]    = (unsigned int) *step_data; step_data++;// get the step
		head_float[0]  =                *time_data; time_data++;// get the time
		for(i=0; i<9; i++)
		{
			head_float[i+1] = (*box_data) / scaling;
			box_data++;
		}
		xdrfile_write_int( head_int, 3, xfp);
		xdrfile_write_float( head_float, 10, xfp);
		if(verb)
			printf("\r\t... frame %d of %d", T+1, T_frames);



		///######################## Part copied from xdrfile.c and adapted to fit to numpy ###############################////
		//*

		precision = start_prec;

		//check if file pointer is valid
		if(xfp==NULL)
		{
			PyErr_SetString(PyExc_RuntimeError, "Could not create File.\n");
			return Py_BuildValue("i", -1);
		}
	    bitsizeint[0] = 0;
	    bitsizeint[1] = 0;
	    bitsizeint[2] = 0;

		if(xdrfile_write_int(&N_atoms,1,xfp)==0)
			return -1; // return if we could not write size //
		// Dont bother with compression for three atoms or less //
		if(N_atoms<=9)
	    {
			return xdrfile_write_float(ptr,size3,xfp)/3;
			// return number of coords, not floats //
		}
		// Compression-time if we got here. Write precision first //
		if (precision <= 0)
			precision = 1000;
		xdrfile_write_float(&precision,1,xfp);
		// avoid repeated pointer dereferencing. //
		buf1=xfp->buf1;
		buf2=xfp->buf2;
		// buf2[0-2] are special and do not contain actual data //
		buf2[0] = buf2[1] = buf2[2] = 0;
		minint[0] = minint[1] = minint[2] = INT_MAX;
		maxint[0] = maxint[1] = maxint[2] = INT_MIN;
		prevrun = -1;
		lfp = ptr;
		lip = buf1;
		mindiff = INT_MAX;
		oldlint1 = oldlint2 = oldlint3 = 0;

		//while(lfp < ptr + size3 )
		ptr =PyArray_GETPTR3(x, T, 0, 0);
		precision /= scaling;
		for(N=0; N<N_atoms; N++)
	    {
			// find nearest integer //
			//// X
			lf = ptr[N*X_Nstride]; //[N*X_Nstride + 0*X_Cstride]
			if (lf >= 0.0)
				lf = lf * precision + 0.5;
			else
				lf = lf * precision - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				// scaling would cause overflow //
				PyErr_SetString(PyExc_RuntimeError, "Internal overflow compressing coordinates.\n");
				return Py_BuildValue("i", -1);
			}
			lint1 = lf;
			if (lint1 < minint[0]) minint[0] = lint1;
			if (lint1 > maxint[0]) maxint[0] = lint1;
			*lip++ = lint1;

			//// Y
			lf = ptr[N*X_Nstride + X_Cstride]; //[N*X_Nstride + 1*X_Cstride]
			if (lf >= 0.0)
				lf = lf * precision + 0.5;
			else
				lf = lf * precision - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				// scaling would cause overflow */
				PyErr_SetString(PyExc_RuntimeError, "Internal overflow compressing coordinates.\n");
				return Py_BuildValue("i", -1);
			}
			lint2 = lf;
			if (lint2 < minint[1]) minint[1] = lint2;
			if (lint2 > maxint[1]) maxint[1] = lint2;
			*lip++ = lint2;

			//// Z
			lf = ptr[N*X_Nstride + 2*X_Cstride];
			if (lf >= 0.0)
				lf = lf * precision + 0.5;
			else
				lf = lf * precision - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				errval=0;
			}
			lint3 = lf;
			if (lint3 < minint[2]) minint[2] = lint3;
			if (lint3 > maxint[2]) maxint[2] = lint3;
			*lip++ = lint3;

			// wrap up
			diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
			if (diff < mindiff && N > 3  )// && lfp > ptr + 3)
				mindiff = diff;
			oldlint1 = lint1;
			oldlint2 = lint2;
			oldlint3 = lint3;
		}
		xdrfile_write_int(minint,3,xfp);
		xdrfile_write_int(maxint,3,xfp);

		if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
			(float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
			(float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
			// turning value in unsigned by subtracting minint
			// would cause overflow

			PyErr_SetString(PyExc_RuntimeError, "Internal overflow compressing coordinates.\n");
			return Py_BuildValue("i", -1);
		}
		sizeint[0] = maxint[0] - minint[0]+1;
		sizeint[1] = maxint[1] - minint[1]+1;
		sizeint[2] = maxint[2] - minint[2]+1;

		// check if one of the sizes is to big to be multiplied //
		if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
	    {
			bitsizeint[0] = sizeofint(sizeint[0]);
			bitsizeint[1] = sizeofint(sizeint[1]);
			bitsizeint[2] = sizeofint(sizeint[2]);
			bitsize = 0; // flag the use of large sizes //
		}
	    else
	    {
			bitsize = sizeofints(3, sizeint);
		}
		lip = buf1;
		luip = (unsigned int *) buf1;
		smallidx = FIRSTIDX;
		while (smallidx < LASTIDX && magicints[smallidx] < mindiff)
	    {
			smallidx++;
		}
		xdrfile_write_int(&smallidx,1,xfp);
		tmp=smallidx+8;
		maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
		minidx = maxidx - 8; // often this equal smallidx //
		tmp=smallidx-1;
		tmp= (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
		smaller = magicints[tmp] / 2;
		smallnum = magicints[smallidx] / 2;
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
		larger = magicints[maxidx] / 2;
		i = 0;
		while (i < N_atoms)
	    {
			is_small = 0;
			thiscoord = (int *)(luip) + i * 3;
			if (smallidx < maxidx && i >= 1 &&
				abs(thiscoord[0] - prevcoord[0]) < larger &&
				abs(thiscoord[1] - prevcoord[1]) < larger &&
				abs(thiscoord[2] - prevcoord[2]) < larger) {
				is_smaller = 1;
			}
	        else if (smallidx > minidx)
	        {
				is_smaller = -1;
			}
	        else
	        {
				is_smaller = 0;
			}
			if (i + 1 < N_atoms)
	        {
				if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
					abs(thiscoord[1] - thiscoord[4]) < smallnum &&
					abs(thiscoord[2] - thiscoord[5]) < smallnum)
	            {
					// interchange first with second atom for better
					// compression of water molecules

					tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
					thiscoord[3] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
					thiscoord[4] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
					thiscoord[5] = tmp;
					is_small = 1;
				}
			}
			tmpcoord[0] = thiscoord[0] - minint[0];
			tmpcoord[1] = thiscoord[1] - minint[1];
			tmpcoord[2] = thiscoord[2] - minint[2];
			if (bitsize == 0)
	        {
				encodebits(buf2, bitsizeint[0], tmpcoord[0]);
				encodebits(buf2, bitsizeint[1], tmpcoord[1]);
				encodebits(buf2, bitsizeint[2], tmpcoord[2]);
			}
	        else
	        {
				encodeints(buf2, 3, bitsize, sizeint, tmpcoord);
			}
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
			thiscoord = thiscoord + 3;
			i++;

			run = 0;
			if (is_small == 0 && is_smaller == -1)
				is_smaller = 0;
			while (is_small && run < 8*3)
	        {
				tmpsum=0;
				for(j=0;j<3;j++)
	            {
					tmp=thiscoord[j] - prevcoord[j];
					tmpsum+=tmp*tmp;
				}
				if (is_smaller == -1 && tmpsum >= smaller * smaller)
	            {
					is_smaller = 0;
				}

				tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
				tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
				tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;

				prevcoord[0] = thiscoord[0];
				prevcoord[1] = thiscoord[1];
				prevcoord[2] = thiscoord[2];

				i++;
				thiscoord = thiscoord + 3;
				is_small = 0;
				if (i < size &&
					abs(thiscoord[0] - prevcoord[0]) < smallnum &&
					abs(thiscoord[1] - prevcoord[1]) < smallnum &&
					abs(thiscoord[2] - prevcoord[2]) < smallnum)
	            {
					is_small = 1;
				}
			}
			if (run != prevrun || is_smaller != 0)
	        {
				prevrun = run;
				encodebits(buf2, 1, 1); // flag the change in run-length //
				encodebits(buf2, 5, run+is_smaller+1);
			}
	        else
	        {
				encodebits(buf2, 1, 0); // flag the fact that runlength did not change //
			}
			for (k=0; k < run; k+=3)
	        {
				encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]);
			}
			if (is_smaller != 0)
	        {
				smallidx += is_smaller;
				if (is_smaller < 0)
	            {
					smallnum = smaller;
					smaller = magicints[smallidx-1] / 2;
				}
	            else
	            {
					smaller = smallnum;
					smallnum = magicints[smallidx] / 2;
				}
				sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
			}
		}
		if (buf2[1] != 0) buf2[0]++;
		xdrfile_write_int(buf2,1,xfp); // buf2[0] holds the length in bytes //
		tmp=xdrfile_write_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp);
		if(!tmp==(unsigned int)buf2[0])
		{
			PyErr_SetString(PyExc_RuntimeError, "Error writing to XTC file.\n");
			return Py_BuildValue("i", -1);
		}

		////// xdrfile.c end
		//*/

	}//end while( T < T_frames )
	if(verb)
		printf(" ... done\n");
	return Py_BuildValue("i", N_atoms);
}






/*
static PyObject *xdrfile_writeI(PyObject *self, PyObject *args)
{
	const char 		*file_path;						//file path of XTC file
	XDRFILE 		*xfp;							//file handle
	XDR_HEAD 		*header;						//file header (with all important info)
	PyArrayObject	*x, *time, *box, *step;				//Python arrays
	float 			*x_data, *time_data, *box_data;	//pointers to tha data of the Python array
	int				X_Tstride, X_Nstride, X_Cstride; //the strides to move through the array
	long int		*step_data;
	float 			scaling=10.0;
	int				precision=1000, start, stop, skip, verb;

	if (!PyArg_ParseTuple(args, "sO!O!O!O!|fii:xdrfile_write",
    	&file_path,
        &PyArray_Type, &x,
        &PyArray_Type, &box,
        &PyArray_Type, &time,
        &PyArray_Type, &step,
        &scaling, &verb, precision)			)  return Py_BuildValue("i", 0);


	unsigned int	   T_frames, N_atoms;
	long unsigned int  T=0,  N=0;

	T_frames 	= (int) PyArray_DIM(x, 0); //number of time frames
	N_atoms		= (int) PyArray_DIM(x, 1); //number of atoms
	time_data	= PyArray_DATA( time );
	box_data    = PyArray_DATA( box );
	step_data   = PyArray_DATA( step );
	X_Tstride   = PyArray_STRIDE(x, 0); //stride along time frame axis
	X_Nstride   = PyArray_STRIDE(x, 1); //       along atoms axis
	X_Cstride   = PyArray_STRIDE(x, 2); //       along coordinates axis (X,Y,Z)

	printf("Contains %d Frames with %d Atoms \n", T_frames, N_atoms);

	xfp = xdrfile_open(file_path, "w+"); //open file to write...
	/// Part copied from xdrfile.c and adapted to fit to numpy ////

	//definitions
	int minint[3], maxint[3], mindiff, *lip, diff;
	int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
	int minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size, size3, *luip;
	int k, *buf1, *buf2;
	int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
	float *ptr, *lfp, lf;
	int tmp, tmpsum, *thiscoord,  prevcoord[3];
	unsigned int tmpcoord[30];
	int errval=1;
	unsigned int bitsize;

	//for the headder
	int   head_int[3];
	float head_float[10];
	head_int[0] = 1995;
	head_int[1] = N_atoms;


	// definitions that need to be done just once
	size3=3*N_atoms;
	if(size3>xfp->buf1size)
    {
		if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL)
        {
			PyErr_SetString(PyExc_BufferError, "Cannot allocate memory for compressing coordinates.\n");
			return Py_BuildValue("i", -1);
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL)
        {
			PyErr_SetString(PyExc_BufferError, "Cannot allocate memory for compressing coordinates.\n");
			return Py_BuildValue("i", -1);
		}
	}


	while( T < T_frames )
	{

		//Write the header
		head_int[2]    = (unsigned int) *step_data; step_data++;// get the step
		head_float[0]  =                *time_data; time_data++;// get the time
		for(i=0; i<9; i++)
		{
			head_float[i+1] = *box_data / scaling;
			box_data++;
		}
		xdrfile_write_int( head_int, 3, xfp);
		xdrfile_write_float( head_float, 10, xfp);
		if(verb)
			printf("writing Frame %d of %d", T+1, T_frames);

		if(xfp==NULL)
		{
			PyErr_SetString(PyExc_RuntimeError, "Could not create File.\n");
			return Py_BuildValue("i", -1);
		}

		ptr  = PyArray_GETPTR1(x, T);
		size = PyArray_DIM(x, 1);
		T++;


		printf("  first floats: %f, %f, %f, %f\n", *ptr, *(ptr+1), *(ptr+2), *(ptr+3) );


		/// Code from xdrfile.c ///////////////////
	    bitsizeint[0] = 0;
	    bitsizeint[1] = 0;
	    bitsizeint[2] = 0;

		/*if(size3>xfp->buf1size)
	    {
			if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL)
	        {
				fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
				return -1;
			}
			xfp->buf1size=size3;
			xfp->buf2size=size3*1.2;
			if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL)
	        {
				fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
				return -1;
			}
		} //Todo reinstate this with the right strides...
		if(xdrfile_write_int(&size,1,xfp)==0)
			return -1; /* return if we could not write size
		 Dont bother with compression for three atoms or less
		if(size<=9)
	    {
			return xdrfile_write_float(ptr,size3,xfp)/3;
			// return number of coords, not floats
		}
		//Compression-time if we got here. Write precision first
		if (precision <= 0)
			precision = 1000;
		xdrfile_write_float(&precision,1,xfp);
		// avoid repeated pointer dereferencing.
		buf1=xfp->buf1;
		buf2=xfp->buf2;
		// buf2[0-2] are special and do not contain actual data
		buf2[0] = buf2[1] = buf2[2] = 0;
		minint[0] = minint[1] = minint[2] = INT_MAX;
		maxint[0] = maxint[1] = maxint[2] = INT_MIN;
		prevrun = -1;
		lfp = ptr;
		lip = buf1;
		mindiff = INT_MAX;
		oldlint1 = oldlint2 = oldlint3 = 0;
		while(lfp < ptr + X_Tstride*T + X_Nstride*(size3/3) + X_Cstride*(size3%3) )
	    {
			// find nearest integer
			if (*lfp >= 0.0)
				lf = *lfp  + 0.5;
			else
				lf = *lfp  - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				// scaling would cause overflow
				fprintf(stderr,"Internal overflow compressing coordinates.\n");
				errval=0;
			}
			lint1 = lf;
			if (lint1 < minint[0]) minint[0] = lint1;
			if (lint1 > maxint[0]) maxint[0] = lint1;
			*lip++ = lint1;
			lfp = ptr + X_Tstride * T + X_Nstride*(N/3) + X_Cstride*(N%3);
			N++;
			if (*lfp >= 0.0)
				lf = *lfp * precision + 0.5;
			else
				lf = *lfp * precision - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				// scaling would cause overflow
				fprintf(stderr,"Internal overflow compressing coordinates.\n");
				errval=0;
			}
			lint2 = lf;
			if (lint2 < minint[1]) minint[1] = lint2;
			if (lint2 > maxint[1]) maxint[1] = lint2;
			*lip++ = lint2;
			lfp = ptr + X_Tstride*T + X_Nstride*(N/3) + X_Cstride*(N%3);
			N++;
			if (*lfp >= 0.0)
				lf = *lfp * precision + 0.5;
			else
				lf = *lfp * precision - 0.5;
			if (fabs(lf) > INT_MAX-2)
	        {
				errval=0;
			}
			lint3 = lf;
			if (lint3 < minint[2]) minint[2] = lint3;
			if (lint3 > maxint[2]) maxint[2] = lint3;
			*lip++ = lint3;
			lfp = ptr + X_Tstride*T + X_Nstride*(N/3) + X_Cstride*(N%3);
			N++;
			diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
			if (diff < mindiff && lfp > ptr + 3)
				mindiff = diff;
			oldlint1 = lint1;
			oldlint2 = lint2;
			oldlint3 = lint3;
		}
		xdrfile_write_int(minint,3,xfp);
		xdrfile_write_int(maxint,3,xfp);

		if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
			(float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
			(float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
			// turning value in unsigned by subtracting minint
			// would cause overflow

			fprintf(stderr,"Internal overflow compressing coordinates.\n");
			errval=0;
		}
		sizeint[0] = maxint[0] - minint[0]+1;
		sizeint[1] = maxint[1] - minint[1]+1;
		sizeint[2] = maxint[2] - minint[2]+1;

		// check if one of the sizes is to big to be multiplied
		if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
	    {
			bitsizeint[0] = sizeofint(sizeint[0]);
			bitsizeint[1] = sizeofint(sizeint[1]);
			bitsizeint[2] = sizeofint(sizeint[2]);
			bitsize = 0; // flag the use of large sizes
		}
	    else
	    {
			bitsize = sizeofints(3, sizeint);
		}
		lip = buf1;
		luip = (unsigned int *) buf1;
		smallidx = FIRSTIDX;
		while (smallidx < LASTIDX && magicints[smallidx] < mindiff)
	    {
			smallidx++;
		}
		xdrfile_write_int(&smallidx,1,xfp);
		tmp=smallidx+8;
		maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
		minidx = maxidx - 8; // often this equal smallidx
		tmp=smallidx-1;
		tmp= (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
		smaller = magicints[tmp] / 2;
		smallnum = magicints[smallidx] / 2;
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
		larger = magicints[maxidx] / 2;
		i = 0;
		while (i < size)
	    {
			is_small = 0;
			thiscoord = (int *)(luip) + i * 3;
			if (smallidx < maxidx && i >= 1 &&
				abs(thiscoord[0] - prevcoord[0]) < larger &&
				abs(thiscoord[1] - prevcoord[1]) < larger &&
				abs(thiscoord[2] - prevcoord[2]) < larger) {
				is_smaller = 1;
			}
	        else if (smallidx > minidx)
	        {
				is_smaller = -1;
			}
	        else
	        {
				is_smaller = 0;
			}
			if (i + 1 < size)
	        {
				if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
					abs(thiscoord[1] - thiscoord[4]) < smallnum &&
					abs(thiscoord[2] - thiscoord[5]) < smallnum)
	            {
    				// interchange first with second atom for better
					// compression of water molecules

					tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
					thiscoord[3] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
					thiscoord[4] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
					thiscoord[5] = tmp;
					is_small = 1;
				}
			}
			tmpcoord[0] = thiscoord[0] - minint[0];
			tmpcoord[1] = thiscoord[1] - minint[1];
			tmpcoord[2] = thiscoord[2] - minint[2];
			if (bitsize == 0)
	        {
				encodebits(buf2, bitsizeint[0], tmpcoord[0]);
				encodebits(buf2, bitsizeint[1], tmpcoord[1]);
				encodebits(buf2, bitsizeint[2], tmpcoord[2]);
			}
	        else
	        {
				encodeints(buf2, 3, bitsize, sizeint, tmpcoord);
			}
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
			thiscoord = thiscoord + 3;
			i++;

			run = 0;
			if (is_small == 0 && is_smaller == -1)
				is_smaller = 0;
			while (is_small && run < 8*3)
	        {
				tmpsum=0;
				for(j=0;j<3;j++)
	            {
					tmp=thiscoord[j] - prevcoord[j];
					tmpsum+=tmp*tmp;
				}
				if (is_smaller == -1 && tmpsum >= smaller * smaller)
	            {
					is_smaller = 0;
				}

				tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
				tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
				tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;

				prevcoord[0] = thiscoord[0];
				prevcoord[1] = thiscoord[1];
				prevcoord[2] = thiscoord[2];

				i++;
				thiscoord = thiscoord + 3;
				is_small = 0;
				if (i < size &&
					abs(thiscoord[0] - prevcoord[0]) < smallnum &&
					abs(thiscoord[1] - prevcoord[1]) < smallnum &&
					abs(thiscoord[2] - prevcoord[2]) < smallnum)
	            {
					is_small = 1;
				}
			}
			if (run != prevrun || is_smaller != 0)
	        {
				prevrun = run;
				encodebits(buf2, 1, 1); // flag the change in run-length
				encodebits(buf2, 5, run+is_smaller+1);
			}
	        else
	        {
				encodebits(buf2, 1, 0); // flag the fact that runlength did not change
			}
			for (k=0; k < run; k+=3)
	        {
				encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]);
			}
			if (is_smaller != 0)
	        {
				smallidx += is_smaller;
				if (is_smaller < 0)
	            {
					smallnum = smaller;
					smaller = magicints[smallidx-1] / 2;
				}
	            else
	            {
					smaller = smallnum;
					smallnum = magicints[smallidx] / 2;
				}
				sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
			}
		}
		if (buf2[1] != 0) buf2[0]++;
		xdrfile_write_int(buf2,1,xfp); // buf2[0] holds the length in bytes
		tmp=xdrfile_write_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp);
		if(!tmp==(unsigned int)buf2[0])
		{
			PyErr_SetString(PyExc_Warning, "Internal error while writing coordinates.\n");
			return Py_BuildValue("i", -1);
		}
	}//end while( T < T_frames )
	if(verb)
		printf("...done\n");
	return Py_BuildValue("i", N_atoms);
}*/


////////////////////////OLD TEST STUFF ////////////////////////
///////////////////////////////////////////////////////////////
static PyObject *xdrfile_readI(PyObject *self, PyObject *args)
{
//self.filename, DatObj.x, DatObj.time, DatObj.BOX, start, stop, skip, verbose=True

	const char *file_path;				//file path of XTC file
	long file_byte_len, temp_file_pos;	//the length of the file in bytes
	XDRFILE* xfp;						//file handle
	XDR_HEAD* header;					//file header (with all important info)

	PyArrayObject *x, *time, *box, *ndx;//Python arrays
	float *x_data, *time_data, *box_data;//poniter to the data of the python arrays
	int *ndx_data, start, stop, skip, verb;		//where to start end and how many fraes to skip (and if verbose or not)

	//decoding stuff
	//int sizeint[3];					//max size for integers of each coordiante
	int T=0, Tr=0;						//Frame counter looping variable; Read frame to adress the arrays
	//int First_run=0;					//indicates that this is the first time in the loop!

	//// Standard lib initializations
	int *minint, *maxint, *lip;
	int smallidx, minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int k, *buf1, *buf2, lsize, flag;
	int smallnum, smaller, larger, i, is_smaller, run;
	float *lfp, inv_precision;
	int tmp, *thiscoord,  prevcoord[3];
	unsigned int bitsize;
	////////////////

	int  N_stride, ndx_stride;
	//float precision;



	//buffer
	int *buf, buf_size=0;


	if (!PyArg_ParseTuple(args, "sO!O!O!O!iiii",
    	&file_path,
        &PyArray_Type, &x,
        &PyArray_Type, &box,
        &PyArray_Type, &time,
        &PyArray_Type, &ndx,
        &start, &stop, &skip,
        &verb					))  return Py_BuildValue("i", 0);

	//opening file handle...
	xfp = xdrfile_open(file_path, "r");
	//...and getting length of the file
	fseek(xfp->fp, 0L, SEEK_END);
	file_byte_len = ftell(xfp->fp);
	fseek(xfp->fp, 0L, SEEK_SET);

	//making space for the header
	header = (XDR_HEAD*) malloc( sizeof(XDR_HEAD) );


	//getting access to the Data of the arrays
	//x_data    = x->data;
	time_data = time->data;
	box_data  = box->data;
	ndx_data  = ndx->data; //we make sure in python that it is sorted!

	N_stride 	= PyArray_STRIDE(x, 0) / 4; // not sure why we divide by 4, but that way it works...



	///Initializations independent of time frame or Atom number (we assume that the number of coordinates won't change with the time frame!)
/*	if( !read_header(xfp, header, file_byte_len) )
		return NULL;				//read the header once to initialize variables
	fseek(xfp->fp, 0L, SEEK_SET);	//go back to the start..

	lsize = header->N1;							//Number of Atoms
	size3 = lsize * 3;							//Number of Floats needed (#Atoms * 3dim)
	inv_precision = 1.0 / header->precision;

	//Reserve Buffer to decompress
		///Buffer1 for decoded Data
	//if(   (  xfp->buf1 = (int *) malloc( sizeof(int)*size3 )  ) == NULL   )
	        {
				PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory for decompressing coordinates.\n");
				return NULL;
			}
	buf1=xfp->buf1;

	//Buffer2 for encoded data
	if(  (xfp->buf2=(int *)malloc(sizeof(int)*size3*1.2)) == NULL  )
	        {
				PyErr_SetString(PyExc_MemoryError, "Cannot allocate memory for decompressing coordinates.\n");
				return -1;
			}
	buf2=xfp->buf2; */




	while(T < start) //skip until we are at the frame we want to start reading from to the start position
	{
		//skip over data
		if( !read_header(xfp, header, file_byte_len) )
			return NULL;
		if( !skip_frame(xfp, header, file_byte_len) )
			return NULL;
		T++;
	}


	while( T < stop ) //now read every go on with the looping but read in every "skip"th frame
	{

		//		Read HEADER
		if( !read_header(xfp, header, file_byte_len) )
		{
			xdrfile_close(xfp); //clean exit ...
			free(header);
			//PyErr_SetString(PyExc_MemoryError, "could not read header");
			return NULL;
		}

		//		TEST FOR SKIP AND DO SO
		if( (T-start) % skip) //if expr==0 --> False
		{
			//skip over data
			if( !skip_frame(xfp, header, file_byte_len) )
			{
				PyErr_SetString(PyExc_MemoryError, "could not jump");
				xdrfile_close(xfp);
				free(header);
				return NULL;
			}
			T++;
			continue; //avoid the else clause...
		}


		// IF WE ARE HERE IT'S READING TIME
		printf("reading frame %i\r",T);

		// write box, and time info to the array
		time_data[Tr] = header->time;
		for(i=0; i<9; i++)
			box_data[ Tr*9 + i ] = header->box[i];


		////////////////////////
		////////////////////////

	    bitsizeint[0] = 0;
	    bitsizeint[1] = 0;
	    bitsizeint[2] = 0;

		//if(xfp==NULL || ptr==NULL)
		//	return -1;
		//tmp=xdrfile_read_int(&lsize,1,xfp);

		//if(tmp==0)
		//	return -1; /* return if we could not read size */
		//if (*size < lsize)
	    //{
		//	fprintf(stderr, "Requested to decompress %d coords, file contains %d\n",
		//			*size, lsize);
		//	return -1;
		//}


		lsize = header->N1;
		size3 = lsize * 3;
		if(size3>xfp->buf1size)
	    {
			if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL)
	        {
				fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
				return -1;
			}
			xfp->buf1size=size3;
			xfp->buf2size=size3*1.2;
			if(  (xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size)) == NULL  )
	        {
				fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
				return -1;
			}
		}
		/* Dont bother with compression for three atoms or less /
		if(*size<=9)
	    {
			return xdrfile_read_float(ptr,size3,xfp)/3;
			return number of coords, not floats
		}*/
		/* Compression-time if we got here. Read precision first */
		//xdrfile_read_float(precision,1,xfp);

		/* avoid repeated pointer dereferencing. */
		buf1=xfp->buf1;
		buf2=xfp->buf2;
		/* buf2[0-2] are special and do not contain actual data */
		buf2[0] = buf2[1] = buf2[2] = 0;


		//xdrfile_read_int(minint,3,xfp);
		//xdrfile_read_int(maxint,3,xfp);
		maxint = header->maxint;
		minint = header->minint;
		sizeint[0] = maxint[0] - minint[0]+1;
		sizeint[1] = maxint[1] - minint[1]+1;
		sizeint[2] = maxint[2] - minint[2]+1;

		/* check if one of the sizes is to big to be multiplied */
		if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
	    {
			bitsizeint[0] = sizeofint(sizeint[0]);
			bitsizeint[1] = sizeofint(sizeint[1]);
			bitsizeint[2] = sizeofint(sizeint[2]);
			bitsize = 0; /* flag the use of large sizes */
		}
	    else
	    {
			bitsize = sizeofints(3, sizeint);
		}

		//if (xdrfile_read_int(&smallidx,1,xfp) == 0)
		//	return 0; /* not sure what has happened here or why we return... */
		smallidx = header->smallidx;

		tmp=smallidx+8;
		maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
		minidx = maxidx - 8; /* often this equal smallidx */
		tmp = smallidx-1;
		tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
		smaller = magicints[tmp] / 2;
		smallnum = magicints[smallidx] / 2;
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
		larger = magicints[maxidx];

		/* buf2[0] holds the length in bytes */

		//if (xdrfile_read_int(buf2,1,xfp) == 0)
		//	return 0;
		buf2[0] = header->data_len;
		if (xdrfile_read_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp) == 0)
			return 0;
		buf2[0] = buf2[1] = buf2[2] = 0;

		lfp = PyArray_GETPTR2(x, 0, Tr);
		inv_precision = 10.0 / header->precision;
		run = 0;
		i = 0;
		lip = buf1;

		//printf("%i:%i '%i' %f\n", buf2[0], buf2[1], buf2[3], inv_precision);
		while ( i < lsize )
	    {

			//if (Tr == 0)
			//	printf("i: %i --> %i:%i:%u\t",i, buf2[0], buf2[1], buf2[2]);

			thiscoord = (int *)(lip) + i * 3;

			if (bitsize == 0)
	        {
				thiscoord[0] = decodebits(buf2, bitsizeint[0]);
				thiscoord[1] = decodebits(buf2, bitsizeint[1]);
				thiscoord[2] = decodebits(buf2, bitsizeint[2]);
			}
	        else
	        {
				decodeints(buf2, 3, bitsize, sizeint, thiscoord);
			}

			//if (Tr == 0)
			//	printf("after read %i:%i\t", buf2[0], buf2[1]);

			i++;
			thiscoord[0] += minint[0];
			thiscoord[1] += minint[1];
			thiscoord[2] += minint[2];

			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];


			flag = decodebits(buf2, 1);
			is_smaller = 0;
			if (flag == 1)
	        {
				run = decodebits(buf2, 5);
				is_smaller = run % 3;
				run -= is_smaller;
				is_smaller--;
			}


			//if (Tr == 0)
			//	printf("flag: %i run:%i\n", flag, run);


			if (run > 0)
	        {
				thiscoord += 3;
				for (k = 0; k < run; k+=3)
	            {
					decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
					i++;
					thiscoord[0] += prevcoord[0] - smallnum;
					thiscoord[1] += prevcoord[1] - smallnum;
					thiscoord[2] += prevcoord[2] - smallnum;
					if (k == 0) {
						/* interchange first with second atom for better
						 * compression of water molecules
						 */
						tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
						prevcoord[0] = tmp;
						tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
						prevcoord[1] = tmp;
						tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
						prevcoord[2] = tmp;
						//printf("%i: %f\n", i, thiscoord[0] * inv_precision);
						lfp[0] = prevcoord[0] * inv_precision;
						lfp[1] = prevcoord[1] * inv_precision;
						lfp[2] = prevcoord[2] * inv_precision;
						lfp    += N_stride;
					} else {
						prevcoord[0] = thiscoord[0];
						prevcoord[1] = thiscoord[1];
						prevcoord[2] = thiscoord[2];
					}
					//printf("%i: %f\n", i, thiscoord[0] * inv_precision);
					lfp[0] = thiscoord[0] * inv_precision;
					lfp[1] = thiscoord[1] * inv_precision;
					lfp[2] = thiscoord[2] * inv_precision;
					lfp    += N_stride;
				}
			}
	        else
	        {
				lfp[0] = thiscoord[0] * inv_precision;
				lfp[1] = thiscoord[1] * inv_precision;
				lfp[2] = thiscoord[2] * inv_precision;
				lfp    += N_stride;
			}
			smallidx += is_smaller;
			if (is_smaller < 0)
	        {
				smallnum = smaller;

				if (smallidx > FIRSTIDX)
	            {
					smaller = magicints[smallidx - 1] /2;
				}
	            else
	            {
					smaller = 0;
				}
			}
	        else if (is_smaller > 0)
	        {
				smaller = smallnum;
				smallnum = magicints[smallidx] / 2;
			}
			sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
		}

		////////////////////////
		////////////////////////


		T++; Tr++; //increment the "frame number" and the "read frames number"
	}
	printf("\nClosing....\n");
	xdrfile_close(xfp);
	free(header);
	return Py_BuildValue("i", 1);
};




