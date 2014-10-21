



<!DOCTYPE html>
<html>
<head>
 <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" >
 <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1" >
 
 <meta name="ROBOTS" content="NOARCHIVE">
 
 <link rel="icon" type="image/vnd.microsoft.icon" href="https://ssl.gstatic.com/codesite/ph/images/phosting.ico">
 
 
 <script type="text/javascript">
 
 
 
 
 var codesite_token = "3Cd3YLziNQwHJ6q0INBaXA2gZls:1366032649547";
 
 
 var CS_env = {"token":"3Cd3YLziNQwHJ6q0INBaXA2gZls:1366032649547","projectName":"pmx","domainName":null,"assetHostPath":"https://ssl.gstatic.com/codesite/ph","loggedInUserEmail":"vytautas.gapsys@gmail.com","profileUrl":"/u/110130407061490526737/","assetVersionPath":"https://ssl.gstatic.com/codesite/ph/14689258884487974863","projectHomeUrl":"/p/pmx","relativeBaseUrl":""};
 var _gaq = _gaq || [];
 _gaq.push(
 ['siteTracker._setAccount', 'UA-18071-1'],
 ['siteTracker._trackPageview']);
 
 (function() {
 var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
 ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
 (document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(ga);
 })();
 
 </script>
 
 
 <title>DTI_analysis.py - 
 pmx -
 
 
 python library and tools for computational and structural biophysics - Google Project Hosting
 </title>
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/14689258884487974863/css/core.css">
 
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/14689258884487974863/css/ph_detail.css" >
 
 
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/14689258884487974863/css/d_sb.css" >
 
 
 
<!--[if IE]>
 <link type="text/css" rel="stylesheet" href="https://ssl.gstatic.com/codesite/ph/14689258884487974863/css/d_ie.css" >
<![endif]-->
 <style type="text/css">
 .menuIcon.off { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 -42px }
 .menuIcon.on { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 -28px }
 .menuIcon.down { background: no-repeat url(https://ssl.gstatic.com/codesite/ph/images/dropdown_sprite.gif) 0 0; }
 
 
 
  tr.inline_comment {
 background: #fff;
 vertical-align: top;
 }
 div.draft, div.published {
 padding: .3em;
 border: 1px solid #999; 
 margin-bottom: .1em;
 font-family: arial, sans-serif;
 max-width: 60em;
 }
 div.draft {
 background: #ffa;
 } 
 div.published {
 background: #e5ecf9;
 }
 div.published .body, div.draft .body {
 padding: .5em .1em .1em .1em;
 max-width: 60em;
 white-space: pre-wrap;
 white-space: -moz-pre-wrap;
 white-space: -pre-wrap;
 white-space: -o-pre-wrap;
 word-wrap: break-word;
 font-size: 1em;
 }
 div.draft .actions {
 margin-left: 1em;
 font-size: 90%;
 }
 div.draft form {
 padding: .5em .5em .5em 0;
 }
 div.draft textarea, div.published textarea {
 width: 95%;
 height: 10em;
 font-family: arial, sans-serif;
 margin-bottom: .5em;
 }

 
 .nocursor, .nocursor td, .cursor_hidden, .cursor_hidden td {
 background-color: white;
 height: 2px;
 }
 .cursor, .cursor td {
 background-color: darkblue;
 height: 2px;
 display: '';
 }
 
 
.list {
 border: 1px solid white;
 border-bottom: 0;
}

 
 </style>
</head>
<body class="t4">
<script type="text/javascript">
 window.___gcfg = {lang: 'en'};
 (function() 
 {var po = document.createElement("script");
 po.type = "text/javascript"; po.async = true;po.src = "https://apis.google.com/js/plusone.js";
 var s = document.getElementsByTagName("script")[0];
 s.parentNode.insertBefore(po, s);
 })();
</script>
<div class="headbg">

 <div id="gaia">
 

 <span>
 
 
 
 <b>vytautas.gapsys@gmail.com</b>
 
 
 | <a href="/u/110130407061490526737/" id="projects-dropdown" onclick="return false;"
 ><u>My favorites</u> <small>&#9660;</small></a>
 | <a href="/u/110130407061490526737/" onclick="_CS_click('/gb/ph/profile');"
 title="Profile, Updates, and Settings"
 ><u>Profile</u></a>
 | <a href="https://www.google.com/accounts/Logout?continue=https%3A%2F%2Fcode.google.com%2Fp%2Fpmx%2Fsource%2Fbrowse%2Fscripts%2FDTI_analysis.py" 
 onclick="_CS_click('/gb/ph/signout');"
 ><u>Sign out</u></a>
 
 </span>

 </div>

 <div class="gbh" style="left: 0pt;"></div>
 <div class="gbh" style="right: 0pt;"></div>
 
 
 <div style="height: 1px"></div>
<!--[if lte IE 7]>
<div style="text-align:center;">
Your version of Internet Explorer is not supported. Try a browser that
contributes to open source, such as <a href="http://www.firefox.com">Firefox</a>,
<a href="http://www.google.com/chrome">Google Chrome</a>, or
<a href="http://code.google.com/chrome/chromeframe/">Google Chrome Frame</a>.
</div>
<![endif]-->



 <table style="padding:0px; margin: 0px 0px 10px 0px; width:100%" cellpadding="0" cellspacing="0"
 itemscope itemtype="http://schema.org/CreativeWork">
 <tr style="height: 58px;">
 
 
 
 <td id="plogo">
 <link itemprop="url" href="/p/pmx">
 <a href="/p/pmx/">
 
 
 <img src="/p/pmx/logo?cct=1355339915"
 alt="Logo" itemprop="image">
 
 </a>
 </td>
 
 <td style="padding-left: 0.5em">
 
 <div id="pname">
 <a href="/p/pmx/"><span itemprop="name">pmx</span></a>
 </div>
 
 <div id="psum">
 <a id="project_summary_link"
 href="/p/pmx/"><span itemprop="description">python library and tools for computational and structural biophysics</span></a>
 
 </div>
 
 
 </td>
 <td style="white-space:nowrap;text-align:right; vertical-align:bottom;">
 
 <form action="/hosting/search">
 <input size="30" name="q" value="" type="text">
 
 <input type="submit" name="projectsearch" value="Search projects" >
 </form>
 
 </tr>
 </table>

</div>

 
<div id="mt" class="gtb"> 
 <a href="/p/pmx/" class="tab ">Project&nbsp;Home</a>
 
 
 
 
 <a href="/p/pmx/downloads/list" class="tab ">Downloads</a>
 
 
 
 
 
 <a href="/p/pmx/w/list" class="tab ">Wiki</a>
 
 
 
 
 
 <a href="/p/pmx/issues/list"
 class="tab ">Issues</a>
 
 
 
 
 
 <a href="/p/pmx/source/checkout"
 class="tab active">Source</a>
 
 
 
 
 
 
 
 
 <div class=gtbc></div>
</div>
<table cellspacing="0" cellpadding="0" width="100%" align="center" border="0" class="st">
 <tr>
 
 
 
 
 
 
 <td class="subt">
 <div class="st2">
 <div class="isf">
 
 <form action="/p/pmx/source/browse" style="display: inline">
 
 Repository:
 <select name="repo" id="repo" style="font-size: 92%" onchange="submit()">
 <option value="default">default</option><option value="wiki">wiki</option>
 </select>
 </form>
 
 


 <span class="inst1"><a href="/p/pmx/source/checkout">Checkout</a></span> &nbsp;
 <span class="inst2"><a href="/p/pmx/source/browse/">Browse</a></span> &nbsp;
 <span class="inst3"><a href="/p/pmx/source/list">Changes</a></span> &nbsp;
 <span class="inst4"><a href="/p/pmx/source/clones">Clones</a></span> &nbsp; 
 
 
 
 
 <a href="/p/pmx/issues/entry?show=review&former=sourcelist">Request code review</a>
 
 
 
 </form>
 <script type="text/javascript">
 
 function codesearchQuery(form) {
 var query = document.getElementById('q').value;
 if (query) { form.action += '%20' + query; }
 }
 </script>
 </div>
</div>

 </td>
 
 
 
 <td align="right" valign="top" class="bevel-right"></td>
 </tr>
</table>


<script type="text/javascript">
 var cancelBubble = false;
 function _go(url) { document.location = url; }
</script>
<div id="maincol"
 
>

 




<div class="collapse">
<div id="colcontrol">
<style type="text/css">
 #file_flipper { white-space: nowrap; padding-right: 2em; }
 #file_flipper.hidden { display: none; }
 #file_flipper .pagelink { color: #0000CC; text-decoration: underline; }
 #file_flipper #visiblefiles { padding-left: 0.5em; padding-right: 0.5em; }
</style>
<table id="nav_and_rev" class="list"
 cellpadding="0" cellspacing="0" width="100%">
 <tr>
 
 <td nowrap="nowrap" class="src_crumbs src_nav" width="33%">
 <strong class="src_nav">Source path:&nbsp;</strong>
 <span id="crumb_root">
 
 <a href="/p/pmx/source/browse/">git</a>/&nbsp;</span>
 <span id="crumb_links" class="ifClosed"><a href="/p/pmx/source/browse/scripts/">scripts</a><span class="sp">/&nbsp;</span>DTI_analysis.py</span>
 
 
 
 
 
 <form class="src_nav">
 
 <span class="sourcelabel"><strong>Branch:</strong>
 <select id="branch_select" name="name" onchange="submit()">
 
 <option value="David"
 >
 David
 </option>
 
 <option value="Upgradegrom4.6"
 >
 Upgradegrom4.6
 </option>
 
 <option value="master"
 selected>
 master
 </option>
 
 
 </select>
 </span>
 </form>
 
 
 
 
 



 </td>
 
 
 <td nowrap="nowrap" width="33%" align="center">
 <a href="/p/pmx/source/browse/scripts/DTI_analysis.py?edit=1"
 ><img src="https://ssl.gstatic.com/codesite/ph/images/pencil-y14.png"
 class="edit_icon">Edit file</a>
 </td>
 
 
 <td nowrap="nowrap" width="33%" align="right">
 <table cellpadding="0" cellspacing="0" style="font-size: 100%"><tr>
 
 
 <td class="flipper">
 <ul class="leftside">
 
 <li><a href="/p/pmx/source/browse/scripts/DTI_analysis.py?r=a2102ac8113476c16e34079d2812b130339f54bb" title="Previous">&lsaquo;a2102ac81134</a></li>
 
 </ul>
 </td>
 
 <td class="flipper"><b>82a17baf41be</b></td>
 
 </tr></table>
 </td> 
 </tr>
</table>

<div class="fc">
 
 
 
<style type="text/css">
.undermouse span {
 background-image: url(https://ssl.gstatic.com/codesite/ph/images/comments.gif); }
</style>
<table class="opened" id="review_comment_area"
onmouseout="gutterOut()"><tr>
<td id="nums">
<pre><table width="100%"><tr class="nocursor"><td></td></tr></table></pre>
<pre><table width="100%" id="nums_table_0"><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_1"

 onmouseover="gutterOver(1)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',1);">&nbsp;</span
></td><td id="1"><a href="#1">1</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_2"

 onmouseover="gutterOver(2)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',2);">&nbsp;</span
></td><td id="2"><a href="#2">2</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_3"

 onmouseover="gutterOver(3)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',3);">&nbsp;</span
></td><td id="3"><a href="#3">3</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_4"

 onmouseover="gutterOver(4)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',4);">&nbsp;</span
></td><td id="4"><a href="#4">4</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_5"

 onmouseover="gutterOver(5)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',5);">&nbsp;</span
></td><td id="5"><a href="#5">5</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_6"

 onmouseover="gutterOver(6)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',6);">&nbsp;</span
></td><td id="6"><a href="#6">6</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_7"

 onmouseover="gutterOver(7)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',7);">&nbsp;</span
></td><td id="7"><a href="#7">7</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_8"

 onmouseover="gutterOver(8)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',8);">&nbsp;</span
></td><td id="8"><a href="#8">8</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_9"

 onmouseover="gutterOver(9)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',9);">&nbsp;</span
></td><td id="9"><a href="#9">9</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_10"

 onmouseover="gutterOver(10)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',10);">&nbsp;</span
></td><td id="10"><a href="#10">10</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_11"

 onmouseover="gutterOver(11)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',11);">&nbsp;</span
></td><td id="11"><a href="#11">11</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_12"

 onmouseover="gutterOver(12)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',12);">&nbsp;</span
></td><td id="12"><a href="#12">12</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_13"

 onmouseover="gutterOver(13)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',13);">&nbsp;</span
></td><td id="13"><a href="#13">13</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_14"

 onmouseover="gutterOver(14)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',14);">&nbsp;</span
></td><td id="14"><a href="#14">14</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_15"

 onmouseover="gutterOver(15)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',15);">&nbsp;</span
></td><td id="15"><a href="#15">15</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_16"

 onmouseover="gutterOver(16)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',16);">&nbsp;</span
></td><td id="16"><a href="#16">16</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_17"

 onmouseover="gutterOver(17)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',17);">&nbsp;</span
></td><td id="17"><a href="#17">17</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_18"

 onmouseover="gutterOver(18)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',18);">&nbsp;</span
></td><td id="18"><a href="#18">18</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_19"

 onmouseover="gutterOver(19)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',19);">&nbsp;</span
></td><td id="19"><a href="#19">19</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_20"

 onmouseover="gutterOver(20)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',20);">&nbsp;</span
></td><td id="20"><a href="#20">20</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_21"

 onmouseover="gutterOver(21)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',21);">&nbsp;</span
></td><td id="21"><a href="#21">21</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_22"

 onmouseover="gutterOver(22)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',22);">&nbsp;</span
></td><td id="22"><a href="#22">22</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_23"

 onmouseover="gutterOver(23)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',23);">&nbsp;</span
></td><td id="23"><a href="#23">23</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_24"

 onmouseover="gutterOver(24)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',24);">&nbsp;</span
></td><td id="24"><a href="#24">24</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_25"

 onmouseover="gutterOver(25)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',25);">&nbsp;</span
></td><td id="25"><a href="#25">25</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_26"

 onmouseover="gutterOver(26)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',26);">&nbsp;</span
></td><td id="26"><a href="#26">26</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_27"

 onmouseover="gutterOver(27)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',27);">&nbsp;</span
></td><td id="27"><a href="#27">27</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_28"

 onmouseover="gutterOver(28)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',28);">&nbsp;</span
></td><td id="28"><a href="#28">28</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_29"

 onmouseover="gutterOver(29)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',29);">&nbsp;</span
></td><td id="29"><a href="#29">29</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_30"

 onmouseover="gutterOver(30)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',30);">&nbsp;</span
></td><td id="30"><a href="#30">30</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_31"

 onmouseover="gutterOver(31)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',31);">&nbsp;</span
></td><td id="31"><a href="#31">31</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_32"

 onmouseover="gutterOver(32)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',32);">&nbsp;</span
></td><td id="32"><a href="#32">32</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_33"

 onmouseover="gutterOver(33)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',33);">&nbsp;</span
></td><td id="33"><a href="#33">33</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_34"

 onmouseover="gutterOver(34)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',34);">&nbsp;</span
></td><td id="34"><a href="#34">34</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_35"

 onmouseover="gutterOver(35)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',35);">&nbsp;</span
></td><td id="35"><a href="#35">35</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_36"

 onmouseover="gutterOver(36)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',36);">&nbsp;</span
></td><td id="36"><a href="#36">36</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_37"

 onmouseover="gutterOver(37)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',37);">&nbsp;</span
></td><td id="37"><a href="#37">37</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_38"

 onmouseover="gutterOver(38)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',38);">&nbsp;</span
></td><td id="38"><a href="#38">38</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_39"

 onmouseover="gutterOver(39)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',39);">&nbsp;</span
></td><td id="39"><a href="#39">39</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_40"

 onmouseover="gutterOver(40)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',40);">&nbsp;</span
></td><td id="40"><a href="#40">40</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_41"

 onmouseover="gutterOver(41)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',41);">&nbsp;</span
></td><td id="41"><a href="#41">41</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_42"

 onmouseover="gutterOver(42)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',42);">&nbsp;</span
></td><td id="42"><a href="#42">42</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_43"

 onmouseover="gutterOver(43)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',43);">&nbsp;</span
></td><td id="43"><a href="#43">43</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_44"

 onmouseover="gutterOver(44)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',44);">&nbsp;</span
></td><td id="44"><a href="#44">44</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_45"

 onmouseover="gutterOver(45)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',45);">&nbsp;</span
></td><td id="45"><a href="#45">45</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_46"

 onmouseover="gutterOver(46)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',46);">&nbsp;</span
></td><td id="46"><a href="#46">46</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_47"

 onmouseover="gutterOver(47)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',47);">&nbsp;</span
></td><td id="47"><a href="#47">47</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_48"

 onmouseover="gutterOver(48)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',48);">&nbsp;</span
></td><td id="48"><a href="#48">48</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_49"

 onmouseover="gutterOver(49)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',49);">&nbsp;</span
></td><td id="49"><a href="#49">49</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_50"

 onmouseover="gutterOver(50)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',50);">&nbsp;</span
></td><td id="50"><a href="#50">50</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_51"

 onmouseover="gutterOver(51)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',51);">&nbsp;</span
></td><td id="51"><a href="#51">51</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_52"

 onmouseover="gutterOver(52)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',52);">&nbsp;</span
></td><td id="52"><a href="#52">52</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_53"

 onmouseover="gutterOver(53)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',53);">&nbsp;</span
></td><td id="53"><a href="#53">53</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_54"

 onmouseover="gutterOver(54)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',54);">&nbsp;</span
></td><td id="54"><a href="#54">54</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_55"

 onmouseover="gutterOver(55)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',55);">&nbsp;</span
></td><td id="55"><a href="#55">55</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_56"

 onmouseover="gutterOver(56)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',56);">&nbsp;</span
></td><td id="56"><a href="#56">56</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_57"

 onmouseover="gutterOver(57)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',57);">&nbsp;</span
></td><td id="57"><a href="#57">57</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_58"

 onmouseover="gutterOver(58)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',58);">&nbsp;</span
></td><td id="58"><a href="#58">58</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_59"

 onmouseover="gutterOver(59)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',59);">&nbsp;</span
></td><td id="59"><a href="#59">59</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_60"

 onmouseover="gutterOver(60)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',60);">&nbsp;</span
></td><td id="60"><a href="#60">60</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_61"

 onmouseover="gutterOver(61)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',61);">&nbsp;</span
></td><td id="61"><a href="#61">61</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_62"

 onmouseover="gutterOver(62)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',62);">&nbsp;</span
></td><td id="62"><a href="#62">62</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_63"

 onmouseover="gutterOver(63)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',63);">&nbsp;</span
></td><td id="63"><a href="#63">63</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_64"

 onmouseover="gutterOver(64)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',64);">&nbsp;</span
></td><td id="64"><a href="#64">64</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_65"

 onmouseover="gutterOver(65)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',65);">&nbsp;</span
></td><td id="65"><a href="#65">65</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_66"

 onmouseover="gutterOver(66)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',66);">&nbsp;</span
></td><td id="66"><a href="#66">66</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_67"

 onmouseover="gutterOver(67)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',67);">&nbsp;</span
></td><td id="67"><a href="#67">67</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_68"

 onmouseover="gutterOver(68)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',68);">&nbsp;</span
></td><td id="68"><a href="#68">68</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_69"

 onmouseover="gutterOver(69)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',69);">&nbsp;</span
></td><td id="69"><a href="#69">69</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_70"

 onmouseover="gutterOver(70)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',70);">&nbsp;</span
></td><td id="70"><a href="#70">70</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_71"

 onmouseover="gutterOver(71)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',71);">&nbsp;</span
></td><td id="71"><a href="#71">71</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_72"

 onmouseover="gutterOver(72)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',72);">&nbsp;</span
></td><td id="72"><a href="#72">72</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_73"

 onmouseover="gutterOver(73)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',73);">&nbsp;</span
></td><td id="73"><a href="#73">73</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_74"

 onmouseover="gutterOver(74)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',74);">&nbsp;</span
></td><td id="74"><a href="#74">74</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_75"

 onmouseover="gutterOver(75)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',75);">&nbsp;</span
></td><td id="75"><a href="#75">75</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_76"

 onmouseover="gutterOver(76)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',76);">&nbsp;</span
></td><td id="76"><a href="#76">76</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_77"

 onmouseover="gutterOver(77)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',77);">&nbsp;</span
></td><td id="77"><a href="#77">77</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_78"

 onmouseover="gutterOver(78)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',78);">&nbsp;</span
></td><td id="78"><a href="#78">78</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_79"

 onmouseover="gutterOver(79)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',79);">&nbsp;</span
></td><td id="79"><a href="#79">79</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_80"

 onmouseover="gutterOver(80)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',80);">&nbsp;</span
></td><td id="80"><a href="#80">80</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_81"

 onmouseover="gutterOver(81)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',81);">&nbsp;</span
></td><td id="81"><a href="#81">81</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_82"

 onmouseover="gutterOver(82)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',82);">&nbsp;</span
></td><td id="82"><a href="#82">82</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_83"

 onmouseover="gutterOver(83)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',83);">&nbsp;</span
></td><td id="83"><a href="#83">83</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_84"

 onmouseover="gutterOver(84)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',84);">&nbsp;</span
></td><td id="84"><a href="#84">84</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_85"

 onmouseover="gutterOver(85)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',85);">&nbsp;</span
></td><td id="85"><a href="#85">85</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_86"

 onmouseover="gutterOver(86)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',86);">&nbsp;</span
></td><td id="86"><a href="#86">86</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_87"

 onmouseover="gutterOver(87)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',87);">&nbsp;</span
></td><td id="87"><a href="#87">87</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_88"

 onmouseover="gutterOver(88)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',88);">&nbsp;</span
></td><td id="88"><a href="#88">88</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_89"

 onmouseover="gutterOver(89)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',89);">&nbsp;</span
></td><td id="89"><a href="#89">89</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_90"

 onmouseover="gutterOver(90)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',90);">&nbsp;</span
></td><td id="90"><a href="#90">90</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_91"

 onmouseover="gutterOver(91)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',91);">&nbsp;</span
></td><td id="91"><a href="#91">91</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_92"

 onmouseover="gutterOver(92)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',92);">&nbsp;</span
></td><td id="92"><a href="#92">92</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_93"

 onmouseover="gutterOver(93)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',93);">&nbsp;</span
></td><td id="93"><a href="#93">93</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_94"

 onmouseover="gutterOver(94)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',94);">&nbsp;</span
></td><td id="94"><a href="#94">94</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_95"

 onmouseover="gutterOver(95)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',95);">&nbsp;</span
></td><td id="95"><a href="#95">95</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_96"

 onmouseover="gutterOver(96)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',96);">&nbsp;</span
></td><td id="96"><a href="#96">96</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_97"

 onmouseover="gutterOver(97)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',97);">&nbsp;</span
></td><td id="97"><a href="#97">97</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_98"

 onmouseover="gutterOver(98)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',98);">&nbsp;</span
></td><td id="98"><a href="#98">98</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_99"

 onmouseover="gutterOver(99)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',99);">&nbsp;</span
></td><td id="99"><a href="#99">99</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_100"

 onmouseover="gutterOver(100)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',100);">&nbsp;</span
></td><td id="100"><a href="#100">100</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_101"

 onmouseover="gutterOver(101)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',101);">&nbsp;</span
></td><td id="101"><a href="#101">101</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_102"

 onmouseover="gutterOver(102)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',102);">&nbsp;</span
></td><td id="102"><a href="#102">102</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_103"

 onmouseover="gutterOver(103)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',103);">&nbsp;</span
></td><td id="103"><a href="#103">103</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_104"

 onmouseover="gutterOver(104)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',104);">&nbsp;</span
></td><td id="104"><a href="#104">104</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_105"

 onmouseover="gutterOver(105)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',105);">&nbsp;</span
></td><td id="105"><a href="#105">105</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_106"

 onmouseover="gutterOver(106)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',106);">&nbsp;</span
></td><td id="106"><a href="#106">106</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_107"

 onmouseover="gutterOver(107)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',107);">&nbsp;</span
></td><td id="107"><a href="#107">107</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_108"

 onmouseover="gutterOver(108)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',108);">&nbsp;</span
></td><td id="108"><a href="#108">108</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_109"

 onmouseover="gutterOver(109)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',109);">&nbsp;</span
></td><td id="109"><a href="#109">109</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_110"

 onmouseover="gutterOver(110)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',110);">&nbsp;</span
></td><td id="110"><a href="#110">110</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_111"

 onmouseover="gutterOver(111)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',111);">&nbsp;</span
></td><td id="111"><a href="#111">111</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_112"

 onmouseover="gutterOver(112)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',112);">&nbsp;</span
></td><td id="112"><a href="#112">112</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_113"

 onmouseover="gutterOver(113)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',113);">&nbsp;</span
></td><td id="113"><a href="#113">113</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_114"

 onmouseover="gutterOver(114)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',114);">&nbsp;</span
></td><td id="114"><a href="#114">114</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_115"

 onmouseover="gutterOver(115)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',115);">&nbsp;</span
></td><td id="115"><a href="#115">115</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_116"

 onmouseover="gutterOver(116)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',116);">&nbsp;</span
></td><td id="116"><a href="#116">116</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_117"

 onmouseover="gutterOver(117)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',117);">&nbsp;</span
></td><td id="117"><a href="#117">117</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_118"

 onmouseover="gutterOver(118)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',118);">&nbsp;</span
></td><td id="118"><a href="#118">118</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_119"

 onmouseover="gutterOver(119)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',119);">&nbsp;</span
></td><td id="119"><a href="#119">119</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_120"

 onmouseover="gutterOver(120)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',120);">&nbsp;</span
></td><td id="120"><a href="#120">120</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_121"

 onmouseover="gutterOver(121)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',121);">&nbsp;</span
></td><td id="121"><a href="#121">121</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_122"

 onmouseover="gutterOver(122)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',122);">&nbsp;</span
></td><td id="122"><a href="#122">122</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_123"

 onmouseover="gutterOver(123)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',123);">&nbsp;</span
></td><td id="123"><a href="#123">123</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_124"

 onmouseover="gutterOver(124)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',124);">&nbsp;</span
></td><td id="124"><a href="#124">124</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_125"

 onmouseover="gutterOver(125)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',125);">&nbsp;</span
></td><td id="125"><a href="#125">125</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_126"

 onmouseover="gutterOver(126)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',126);">&nbsp;</span
></td><td id="126"><a href="#126">126</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_127"

 onmouseover="gutterOver(127)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',127);">&nbsp;</span
></td><td id="127"><a href="#127">127</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_128"

 onmouseover="gutterOver(128)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',128);">&nbsp;</span
></td><td id="128"><a href="#128">128</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_129"

 onmouseover="gutterOver(129)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',129);">&nbsp;</span
></td><td id="129"><a href="#129">129</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_130"

 onmouseover="gutterOver(130)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',130);">&nbsp;</span
></td><td id="130"><a href="#130">130</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_131"

 onmouseover="gutterOver(131)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',131);">&nbsp;</span
></td><td id="131"><a href="#131">131</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_132"

 onmouseover="gutterOver(132)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',132);">&nbsp;</span
></td><td id="132"><a href="#132">132</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_133"

 onmouseover="gutterOver(133)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',133);">&nbsp;</span
></td><td id="133"><a href="#133">133</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_134"

 onmouseover="gutterOver(134)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',134);">&nbsp;</span
></td><td id="134"><a href="#134">134</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_135"

 onmouseover="gutterOver(135)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',135);">&nbsp;</span
></td><td id="135"><a href="#135">135</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_136"

 onmouseover="gutterOver(136)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',136);">&nbsp;</span
></td><td id="136"><a href="#136">136</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_137"

 onmouseover="gutterOver(137)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',137);">&nbsp;</span
></td><td id="137"><a href="#137">137</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_138"

 onmouseover="gutterOver(138)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',138);">&nbsp;</span
></td><td id="138"><a href="#138">138</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_139"

 onmouseover="gutterOver(139)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',139);">&nbsp;</span
></td><td id="139"><a href="#139">139</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_140"

 onmouseover="gutterOver(140)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',140);">&nbsp;</span
></td><td id="140"><a href="#140">140</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_141"

 onmouseover="gutterOver(141)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',141);">&nbsp;</span
></td><td id="141"><a href="#141">141</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_142"

 onmouseover="gutterOver(142)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',142);">&nbsp;</span
></td><td id="142"><a href="#142">142</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_143"

 onmouseover="gutterOver(143)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',143);">&nbsp;</span
></td><td id="143"><a href="#143">143</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_144"

 onmouseover="gutterOver(144)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',144);">&nbsp;</span
></td><td id="144"><a href="#144">144</a></td></tr
><tr id="gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_145"

 onmouseover="gutterOver(145)"

><td><span title="Add comment" onclick="codereviews.startEdit('svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b',145);">&nbsp;</span
></td><td id="145"><a href="#145">145</a></td></tr
></table></pre>
<pre><table width="100%"><tr class="nocursor"><td></td></tr></table></pre>
</td>
<td id="lines">
<pre><table width="100%"><tr class="cursor_stop cursor_hidden"><td></td></tr></table></pre>
<pre class="prettyprint lang-py"><table id="src_table_0"><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_1

 onmouseover="gutterOver(1)"

><td class="source"># pmx  Copyright Notice<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_2

 onmouseover="gutterOver(2)"

><td class="source"># ============================<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_3

 onmouseover="gutterOver(3)"

><td class="source">#<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_4

 onmouseover="gutterOver(4)"

><td class="source"># The pmx source code is copyrighted, but you can freely use and<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_5

 onmouseover="gutterOver(5)"

><td class="source"># copy it as long as you don&#39;t change or remove any of the copyright<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_6

 onmouseover="gutterOver(6)"

><td class="source"># notices.<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_7

 onmouseover="gutterOver(7)"

><td class="source">#<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_8

 onmouseover="gutterOver(8)"

><td class="source"># ----------------------------------------------------------------------<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_9

 onmouseover="gutterOver(9)"

><td class="source"># pmx is Copyright (C) 2006-2013 by Daniel Seeliger<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_10

 onmouseover="gutterOver(10)"

><td class="source">#<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_11

 onmouseover="gutterOver(11)"

><td class="source">#                        All Rights Reserved<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_12

 onmouseover="gutterOver(12)"

><td class="source">#<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_13

 onmouseover="gutterOver(13)"

><td class="source"># Permission to use, copy, modify, distribute, and distribute modified<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_14

 onmouseover="gutterOver(14)"

><td class="source"># versions of this software and its documentation for any purpose and<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_15

 onmouseover="gutterOver(15)"

><td class="source"># without fee is hereby granted, provided that the above copyright<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_16

 onmouseover="gutterOver(16)"

><td class="source"># notice appear in all copies and that both the copyright notice and<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_17

 onmouseover="gutterOver(17)"

><td class="source"># this permission notice appear in supporting documentation, and that<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_18

 onmouseover="gutterOver(18)"

><td class="source"># the name of Daniel Seeliger not be used in advertising or publicity<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_19

 onmouseover="gutterOver(19)"

><td class="source"># pertaining to distribution of the software without specific, written<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_20

 onmouseover="gutterOver(20)"

><td class="source"># prior permission.<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_21

 onmouseover="gutterOver(21)"

><td class="source">#<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_22

 onmouseover="gutterOver(22)"

><td class="source"># DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_23

 onmouseover="gutterOver(23)"

><td class="source"># SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_24

 onmouseover="gutterOver(24)"

><td class="source"># FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_25

 onmouseover="gutterOver(25)"

><td class="source"># SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_26

 onmouseover="gutterOver(26)"

><td class="source"># RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_27

 onmouseover="gutterOver(27)"

><td class="source"># CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_28

 onmouseover="gutterOver(28)"

><td class="source"># CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_29

 onmouseover="gutterOver(29)"

><td class="source"># ----------------------------------------------------------------------<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_30

 onmouseover="gutterOver(30)"

><td class="source">import sys, os<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_31

 onmouseover="gutterOver(31)"

><td class="source">from glob import glob<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_32

 onmouseover="gutterOver(32)"

><td class="source">from numpy import *<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_33

 onmouseover="gutterOver(33)"

><td class="source">from pmx import *<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_34

 onmouseover="gutterOver(34)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_35

 onmouseover="gutterOver(35)"

><td class="source">def read_data(fn, b = 0, e = -1):<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_36

 onmouseover="gutterOver(36)"

><td class="source">    if e == -1: e = 9999999999<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_37

 onmouseover="gutterOver(37)"

><td class="source">    l = open(fn).readlines()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_38

 onmouseover="gutterOver(38)"

><td class="source">    data = []<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_39

 onmouseover="gutterOver(39)"

><td class="source">    for line in l:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_40

 onmouseover="gutterOver(40)"

><td class="source">        if line[0] not in [&#39;@&#39;,&#39;#&#39;]:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_41

 onmouseover="gutterOver(41)"

><td class="source">            entr = line.split()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_42

 onmouseover="gutterOver(42)"

><td class="source">            try:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_43

 onmouseover="gutterOver(43)"

><td class="source">                time = float(entr[0])<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_44

 onmouseover="gutterOver(44)"

><td class="source">                if time &gt; b and time &lt; e:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_45

 onmouseover="gutterOver(45)"

><td class="source">                    data.append( float(entr[1] ) )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_46

 onmouseover="gutterOver(46)"

><td class="source">            except:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_47

 onmouseover="gutterOver(47)"

><td class="source">                pass<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_48

 onmouseover="gutterOver(48)"

><td class="source">#    print &gt;&gt;sys.stderr, &#39;Read file:&#39;, fn, &#39; with %d data points&#39; % len(data)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_49

 onmouseover="gutterOver(49)"

><td class="source">    return data<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_50

 onmouseover="gutterOver(50)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_51

 onmouseover="gutterOver(51)"

><td class="source">def datapoint_from_time(time):<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_52

 onmouseover="gutterOver(52)"

><td class="source">    return time*500<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_53

 onmouseover="gutterOver(53)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_54

 onmouseover="gutterOver(54)"

><td class="source">def block_aver( data, block_size = 1000):<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_55

 onmouseover="gutterOver(55)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_56

 onmouseover="gutterOver(56)"

><td class="source">    total_time = len(data) / 500.<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_57

 onmouseover="gutterOver(57)"

><td class="source">    next_time = block_size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_58

 onmouseover="gutterOver(58)"

><td class="source">    results = []<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_59

 onmouseover="gutterOver(59)"

><td class="source">    offset = 0<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_60

 onmouseover="gutterOver(60)"

><td class="source">    while next_time &lt; total_time:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_61

 onmouseover="gutterOver(61)"

><td class="source">        beg = datapoint_from_time(offset)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_62

 onmouseover="gutterOver(62)"

><td class="source">        end = datapoint_from_time(next_time)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_63

 onmouseover="gutterOver(63)"

><td class="source">        res = average( data[beg:end] )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_64

 onmouseover="gutterOver(64)"

><td class="source">        results.append( (str(offset)+&#39;-&#39;+str(next_time), res ) )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_65

 onmouseover="gutterOver(65)"

><td class="source">        offset = next_time<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_66

 onmouseover="gutterOver(66)"

><td class="source">        next_time += block_size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_67

 onmouseover="gutterOver(67)"

><td class="source">    return results<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_68

 onmouseover="gutterOver(68)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_69

 onmouseover="gutterOver(69)"

><td class="source">def convergence( data, block_size = 1000):<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_70

 onmouseover="gutterOver(70)"

><td class="source">    total_time = len(data) / 500.<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_71

 onmouseover="gutterOver(71)"

><td class="source">    next_time = block_size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_72

 onmouseover="gutterOver(72)"

><td class="source">    results = []<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_73

 onmouseover="gutterOver(73)"

><td class="source">    offset = 0<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_74

 onmouseover="gutterOver(74)"

><td class="source">    while next_time &lt; total_time:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_75

 onmouseover="gutterOver(75)"

><td class="source">        beg = datapoint_from_time(offset)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_76

 onmouseover="gutterOver(76)"

><td class="source">        end = datapoint_from_time(next_time)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_77

 onmouseover="gutterOver(77)"

><td class="source">        res = average( data[beg:end] )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_78

 onmouseover="gutterOver(78)"

><td class="source">        results.append( (next_time, res ) )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_79

 onmouseover="gutterOver(79)"

><td class="source">        next_time += block_size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_80

 onmouseover="gutterOver(80)"

><td class="source">    return results<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_81

 onmouseover="gutterOver(81)"

><td class="source">    <br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_82

 onmouseover="gutterOver(82)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_83

 onmouseover="gutterOver(83)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_84

 onmouseover="gutterOver(84)"

><td class="source">help_text = (&#39;Calculate delta G from multiple DTI runs&#39;,)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_85

 onmouseover="gutterOver(85)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_86

 onmouseover="gutterOver(86)"

><td class="source">options = [<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_87

 onmouseover="gutterOver(87)"

><td class="source">        Option( &quot;-b&quot;, &quot;real&quot;, 500, &quot;Start time [ps]&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_88

 onmouseover="gutterOver(88)"

><td class="source">        Option( &quot;-e&quot;, &quot;real&quot;, -1, &quot;End time[ps]&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_89

 onmouseover="gutterOver(89)"

><td class="source">        Option( &quot;-block1&quot;, &quot;int&quot;, 100, &quot;Time[ps] for block average&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_90

 onmouseover="gutterOver(90)"

><td class="source">        Option( &quot;-block2&quot;, &quot;int&quot;, 500, &quot;Time[ps] for block average&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_91

 onmouseover="gutterOver(91)"

><td class="source">#        Option( &quot;-r2&quot;, &quot;rvec&quot;, [1,2,3], &quot;some vector that does wonderful things and returns always segfaults&quot;)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_92

 onmouseover="gutterOver(92)"

><td class="source">        ]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_93

 onmouseover="gutterOver(93)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_94

 onmouseover="gutterOver(94)"

><td class="source">files = [<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_95

 onmouseover="gutterOver(95)"

><td class="source">    FileOption(&quot;-dgdl&quot;, &quot;r&quot;,[&quot;xvg&quot;],&quot;run&quot;, &quot;Input file with dH/dl values&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_96

 onmouseover="gutterOver(96)"

><td class="source">    FileOption(&quot;-o&quot;, &quot;w&quot;,[&quot;txt&quot;],&quot;results.txt&quot;, &quot;Results&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_97

 onmouseover="gutterOver(97)"

><td class="source">    FileOption(&quot;-oc&quot;, &quot;w&quot;,[&quot;txt&quot;],&quot;convergence.txt&quot;, &quot;text file with mutations to insert&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_98

 onmouseover="gutterOver(98)"

><td class="source">    FileOption(&quot;-ob&quot;, &quot;w&quot;,[&quot;txt&quot;],&quot;block.txt&quot;, &quot;files with block averages&quot;),<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_99

 onmouseover="gutterOver(99)"

><td class="source">    <br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_100

 onmouseover="gutterOver(100)"

><td class="source">]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_101

 onmouseover="gutterOver(101)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_102

 onmouseover="gutterOver(102)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_103

 onmouseover="gutterOver(103)"

><td class="source">cmdl = Commandline( sys.argv, options = options,<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_104

 onmouseover="gutterOver(104)"

><td class="source">                    fileoptions = files,<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_105

 onmouseover="gutterOver(105)"

><td class="source">                    program_desc = help_text,<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_106

 onmouseover="gutterOver(106)"

><td class="source">                    check_for_existing_files = False )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_107

 onmouseover="gutterOver(107)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_108

 onmouseover="gutterOver(108)"

><td class="source">dgdl_file = cmdl[&#39;-dgdl&#39;]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_109

 onmouseover="gutterOver(109)"

><td class="source">start_time = cmdl[&#39;-b&#39;]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_110

 onmouseover="gutterOver(110)"

><td class="source">end_time = cmdl[&#39;-e&#39;]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_111

 onmouseover="gutterOver(111)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_112

 onmouseover="gutterOver(112)"

><td class="source">print &#39;DTI_analysis__&gt; Reading: &#39;, dgdl_file<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_113

 onmouseover="gutterOver(113)"

><td class="source">print &#39;DTI_analysis__&gt; Start time = &#39;, start_time, &#39; End time = &#39;, end_time<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_114

 onmouseover="gutterOver(114)"

><td class="source">data = read_data( dgdl_file, b = start_time, e = end_time )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_115

 onmouseover="gutterOver(115)"

><td class="source">av = average(data)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_116

 onmouseover="gutterOver(116)"

><td class="source">st = std(data)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_117

 onmouseover="gutterOver(117)"

><td class="source">size = len(data)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_118

 onmouseover="gutterOver(118)"

><td class="source">print &#39;DTI_analysis__&gt; &lt;dH/dl&gt; = %8.4f&#39;% av, &#39; | #data points = &#39;, size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_119

 onmouseover="gutterOver(119)"

><td class="source">fp = open(cmdl[&#39;-o&#39;],&#39;w&#39;)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_120

 onmouseover="gutterOver(120)"

><td class="source">print &gt;&gt;fp, av, st, size<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_121

 onmouseover="gutterOver(121)"

><td class="source">fp.close()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_122

 onmouseover="gutterOver(122)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_123

 onmouseover="gutterOver(123)"

><td class="source">block1 = cmdl[&#39;-block1&#39;]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_124

 onmouseover="gutterOver(124)"

><td class="source">fn =os.path.splitext(cmdl[&#39;-ob&#39;])[0]+str(block1)+os.path.splitext(cmdl[&#39;-ob&#39;])[1]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_125

 onmouseover="gutterOver(125)"

><td class="source">print &#39;DTI_analysis__&gt; Block averaging 1: &#39;,  block1<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_126

 onmouseover="gutterOver(126)"

><td class="source">res = block_aver( data, block1 )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_127

 onmouseover="gutterOver(127)"

><td class="source">fp = open(fn,&#39;w&#39;)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_128

 onmouseover="gutterOver(128)"

><td class="source">for a, b in res:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_129

 onmouseover="gutterOver(129)"

><td class="source">    print &gt;&gt;fp, a, b<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_130

 onmouseover="gutterOver(130)"

><td class="source">fp.close()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_131

 onmouseover="gutterOver(131)"

><td class="source">block2 = cmdl[&#39;-block2&#39;]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_132

 onmouseover="gutterOver(132)"

><td class="source">fn =os.path.splitext(cmdl[&#39;-ob&#39;])[0]+str(block2)+os.path.splitext(cmdl[&#39;-ob&#39;])[1]<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_133

 onmouseover="gutterOver(133)"

><td class="source">print &#39;DTI_analysis__&gt; Block averaging 2: &#39;,  block2<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_134

 onmouseover="gutterOver(134)"

><td class="source">res = block_aver( data, block2 )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_135

 onmouseover="gutterOver(135)"

><td class="source">fp = open(fn,&#39;w&#39;)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_136

 onmouseover="gutterOver(136)"

><td class="source">for a, b in res:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_137

 onmouseover="gutterOver(137)"

><td class="source">    print &gt;&gt;fp, a, b<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_138

 onmouseover="gutterOver(138)"

><td class="source">fp.close()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_139

 onmouseover="gutterOver(139)"

><td class="source"><br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_140

 onmouseover="gutterOver(140)"

><td class="source">res = convergence( data, 100 )<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_141

 onmouseover="gutterOver(141)"

><td class="source">fp = open(cmdl[&#39;-oc&#39;],&#39;w&#39;)<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_142

 onmouseover="gutterOver(142)"

><td class="source">for t, r in res:<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_143

 onmouseover="gutterOver(143)"

><td class="source">    print &gt;&gt;fp, t, r<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_144

 onmouseover="gutterOver(144)"

><td class="source">fp.close()<br></td></tr
><tr
id=sl_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_145

 onmouseover="gutterOver(145)"

><td class="source"><br></td></tr
></table></pre>
<pre><table width="100%"><tr class="cursor_stop cursor_hidden"><td></td></tr></table></pre>
</td>
</tr></table>

 
<script type="text/javascript">
 var lineNumUnderMouse = -1;
 
 function gutterOver(num) {
 gutterOut();
 var newTR = document.getElementById('gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_' + num);
 if (newTR) {
 newTR.className = 'undermouse';
 }
 lineNumUnderMouse = num;
 }
 function gutterOut() {
 if (lineNumUnderMouse != -1) {
 var oldTR = document.getElementById(
 'gr_svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b_' + lineNumUnderMouse);
 if (oldTR) {
 oldTR.className = '';
 }
 lineNumUnderMouse = -1;
 }
 }
 var numsGenState = {table_base_id: 'nums_table_'};
 var srcGenState = {table_base_id: 'src_table_'};
 var alignerRunning = false;
 var startOver = false;
 function setLineNumberHeights() {
 if (alignerRunning) {
 startOver = true;
 return;
 }
 numsGenState.chunk_id = 0;
 numsGenState.table = document.getElementById('nums_table_0');
 numsGenState.row_num = 0;
 if (!numsGenState.table) {
 return; // Silently exit if no file is present.
 }
 srcGenState.chunk_id = 0;
 srcGenState.table = document.getElementById('src_table_0');
 srcGenState.row_num = 0;
 alignerRunning = true;
 continueToSetLineNumberHeights();
 }
 function rowGenerator(genState) {
 if (genState.row_num < genState.table.rows.length) {
 var currentRow = genState.table.rows[genState.row_num];
 genState.row_num++;
 return currentRow;
 }
 var newTable = document.getElementById(
 genState.table_base_id + (genState.chunk_id + 1));
 if (newTable) {
 genState.chunk_id++;
 genState.row_num = 0;
 genState.table = newTable;
 return genState.table.rows[0];
 }
 return null;
 }
 var MAX_ROWS_PER_PASS = 1000;
 function continueToSetLineNumberHeights() {
 var rowsInThisPass = 0;
 var numRow = 1;
 var srcRow = 1;
 while (numRow && srcRow && rowsInThisPass < MAX_ROWS_PER_PASS) {
 numRow = rowGenerator(numsGenState);
 srcRow = rowGenerator(srcGenState);
 rowsInThisPass++;
 if (numRow && srcRow) {
 if (numRow.offsetHeight != srcRow.offsetHeight) {
 numRow.firstChild.style.height = srcRow.offsetHeight + 'px';
 }
 }
 }
 if (rowsInThisPass >= MAX_ROWS_PER_PASS) {
 setTimeout(continueToSetLineNumberHeights, 10);
 } else {
 alignerRunning = false;
 if (startOver) {
 startOver = false;
 setTimeout(setLineNumberHeights, 500);
 }
 }
 }
 function initLineNumberHeights() {
 // Do 2 complete passes, because there can be races
 // between this code and prettify.
 startOver = true;
 setTimeout(setLineNumberHeights, 250);
 window.onresize = setLineNumberHeights;
 }
 initLineNumberHeights();
</script>

 
 
 <div id="log">
 <div style="text-align:right">
 <a class="ifCollapse" href="#" onclick="_toggleMeta(this); return false">Show details</a>
 <a class="ifExpand" href="#" onclick="_toggleMeta(this); return false">Hide details</a>
 </div>
 <div class="ifExpand">
 
 
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="changelog">
 <p>Change log</p>
 <div>
 <a href="/p/pmx/source/detail?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b">82a17baf41be</a>
 by Daniel Seeliger &lt;seeliged@bibciw20lx.eu.boehringer.com&gt;
 on Mar 22, 2013
 &nbsp; <a href="/p/pmx/source/diff?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;format=side&amp;path=/scripts/DTI_analysis.py&amp;old_path=/scripts/DTI_analysis.py&amp;old=a2102ac8113476c16e34079d2812b130339f54bb">Diff</a>
 </div>
 <pre>changed pdb format (element column)
</pre>
 </div>
 
 
 
 
 
 
 <script type="text/javascript">
 var detail_url = '/p/pmx/source/detail?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b';
 var publish_url = '/p/pmx/source/detail?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b#publish';
 // describe the paths of this revision in javascript.
 var changed_paths = [];
 var changed_urls = [];
 
 changed_paths.push('/pmx/atom.py');
 changed_urls.push('/p/pmx/source/browse/pmx/atom.py?r\x3d82a17baf41be6e3dd11997ac8d7dff4272c3a37b\x26spec\x3dsvn82a17baf41be6e3dd11997ac8d7dff4272c3a37b');
 
 
 changed_paths.push('/pmx/library.py');
 changed_urls.push('/p/pmx/source/browse/pmx/library.py?r\x3d82a17baf41be6e3dd11997ac8d7dff4272c3a37b\x26spec\x3dsvn82a17baf41be6e3dd11997ac8d7dff4272c3a37b');
 
 
 changed_paths.push('/scripts/DTI_analysis.py');
 changed_urls.push('/p/pmx/source/browse/scripts/DTI_analysis.py?r\x3d82a17baf41be6e3dd11997ac8d7dff4272c3a37b\x26spec\x3dsvn82a17baf41be6e3dd11997ac8d7dff4272c3a37b');
 
 var selected_path = '/scripts/DTI_analysis.py';
 
 
 changed_paths.push('/scripts/DTI_analysis2.py');
 changed_urls.push('/p/pmx/source/browse/scripts/DTI_analysis2.py?r\x3d82a17baf41be6e3dd11997ac8d7dff4272c3a37b\x26spec\x3dsvn82a17baf41be6e3dd11997ac8d7dff4272c3a37b');
 
 
 function getCurrentPageIndex() {
 for (var i = 0; i < changed_paths.length; i++) {
 if (selected_path == changed_paths[i]) {
 return i;
 }
 }
 }
 function getNextPage() {
 var i = getCurrentPageIndex();
 if (i < changed_paths.length - 1) {
 return changed_urls[i + 1];
 }
 return null;
 }
 function getPreviousPage() {
 var i = getCurrentPageIndex();
 if (i > 0) {
 return changed_urls[i - 1];
 }
 return null;
 }
 function gotoNextPage() {
 var page = getNextPage();
 if (!page) {
 page = detail_url;
 }
 window.location = page;
 }
 function gotoPreviousPage() {
 var page = getPreviousPage();
 if (!page) {
 page = detail_url;
 }
 window.location = page;
 }
 function gotoDetailPage() {
 window.location = detail_url;
 }
 function gotoPublishPage() {
 window.location = publish_url;
 }
</script>

 
 <style type="text/css">
 #review_nav {
 border-top: 3px solid white;
 padding-top: 6px;
 margin-top: 1em;
 }
 #review_nav td {
 vertical-align: middle;
 }
 #review_nav select {
 margin: .5em 0;
 }
 </style>
 <div id="review_nav">
 <table><tr><td>Go to:&nbsp;</td><td>
 <select name="files_in_rev" onchange="window.location=this.value">
 
 <option value="/p/pmx/source/browse/pmx/atom.py?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b"
 
 >/pmx/atom.py</option>
 
 <option value="/p/pmx/source/browse/pmx/library.py?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b"
 
 >/pmx/library.py</option>
 
 <option value="/p/pmx/source/browse/scripts/DTI_analysis.py?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b"
 selected="selected"
 >/scripts/DTI_analysis.py</option>
 
 <option value="/p/pmx/source/browse/scripts/DTI_analysis2.py?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&amp;spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b"
 
 >/scripts/DTI_analysis2.py</option>
 
 </select>
 </td></tr></table>
 
 
 <div id="review_instr" class="closed">
 <a class="ifOpened" href="/p/pmx/source/detail?r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b&spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b#publish">Publish your comments</a>
 <div class="ifClosed">Double click a line to add a comment</div>
 </div>
 
 </div>
 
 
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="older_bubble">
 <p>Older revisions</p>
 
 
 <div class="closed" style="margin-bottom:3px;" >
 <a class="ifClosed" onclick="return _toggleHidden(this)"><img src="https://ssl.gstatic.com/codesite/ph/images/plus.gif" ></a>
 <a class="ifOpened" onclick="return _toggleHidden(this)"><img src="https://ssl.gstatic.com/codesite/ph/images/minus.gif" ></a>
 <a href="/p/pmx/source/detail?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&r=a2102ac8113476c16e34079d2812b130339f54bb">a2102ac81134</a>
 by Daniel Seeliger &lt;seeliged@bibciw19lx.eu.boehringer.com&gt;
 on Jan 3, 2013
 &nbsp; <a href="/p/pmx/source/diff?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&r=a2102ac8113476c16e34079d2812b130339f54bb&amp;format=side&amp;path=/scripts/DTI_analysis.py&amp;old_path=/scripts/DTI_analysis.py&amp;old=1bc76bef6ad1f114f81454f2a6ed3f365c62f63e">Diff</a>
 <br>
 <pre class="ifOpened">changed MDP class
</pre>
 </div>
 
 <div class="closed" style="margin-bottom:3px;" >
 <a class="ifClosed" onclick="return _toggleHidden(this)"><img src="https://ssl.gstatic.com/codesite/ph/images/plus.gif" ></a>
 <a class="ifOpened" onclick="return _toggleHidden(this)"><img src="https://ssl.gstatic.com/codesite/ph/images/minus.gif" ></a>
 <a href="/p/pmx/source/detail?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&r=1bc76bef6ad1f114f81454f2a6ed3f365c62f63e">1bc76bef6ad1</a>
 by dseelig &lt;seeliger.biosoft@gmail.com&gt;
 on Nov 26, 2012
 &nbsp; <a href="/p/pmx/source/diff?spec=svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b&r=1bc76bef6ad1f114f81454f2a6ed3f365c62f63e&amp;format=side&amp;path=/scripts/DTI_analysis.py&amp;old_path=/scripts/DTI_analysis.py&amp;old=">Diff</a>
 <br>
 <pre class="ifOpened">initial commit
</pre>
 </div>
 
 
 <a href="/p/pmx/source/list?path=/scripts/DTI_analysis.py&r=82a17baf41be6e3dd11997ac8d7dff4272c3a37b">All revisions of this file</a>
 </div>
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 
 <div class="pmeta_bubble_bg" style="border:1px solid white">
 <div class="round4"></div>
 <div class="round2"></div>
 <div class="round1"></div>
 <div class="box-inner">
 <div id="fileinfo_bubble">
 <p>File info</p>
 
 <div>Size: 4898 bytes,
 145 lines</div>
 
 <div><a href="//pmx.googlecode.com/git/scripts/DTI_analysis.py">View raw file</a></div>
 </div>
 
 </div>
 <div class="round1"></div>
 <div class="round2"></div>
 <div class="round4"></div>
 </div>
 </div>
 </div>


</div>

</div>
</div>

<script src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/prettify/prettify.js"></script>
<script type="text/javascript">prettyPrint();</script>


<script src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/source_file_scripts.js"></script>

 <script type="text/javascript" src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/kibbles.js"></script>
 <script type="text/javascript">
 var lastStop = null;
 var initialized = false;
 
 function updateCursor(next, prev) {
 if (prev && prev.element) {
 prev.element.className = 'cursor_stop cursor_hidden';
 }
 if (next && next.element) {
 next.element.className = 'cursor_stop cursor';
 lastStop = next.index;
 }
 }
 
 function pubRevealed(data) {
 updateCursorForCell(data.cellId, 'cursor_stop cursor_hidden');
 if (initialized) {
 reloadCursors();
 }
 }
 
 function draftRevealed(data) {
 updateCursorForCell(data.cellId, 'cursor_stop cursor_hidden');
 if (initialized) {
 reloadCursors();
 }
 }
 
 function draftDestroyed(data) {
 updateCursorForCell(data.cellId, 'nocursor');
 if (initialized) {
 reloadCursors();
 }
 }
 function reloadCursors() {
 kibbles.skipper.reset();
 loadCursors();
 if (lastStop != null) {
 kibbles.skipper.setCurrentStop(lastStop);
 }
 }
 // possibly the simplest way to insert any newly added comments
 // is to update the class of the corresponding cursor row,
 // then refresh the entire list of rows.
 function updateCursorForCell(cellId, className) {
 var cell = document.getElementById(cellId);
 // we have to go two rows back to find the cursor location
 var row = getPreviousElement(cell.parentNode);
 row.className = className;
 }
 // returns the previous element, ignores text nodes.
 function getPreviousElement(e) {
 var element = e.previousSibling;
 if (element.nodeType == 3) {
 element = element.previousSibling;
 }
 if (element && element.tagName) {
 return element;
 }
 }
 function loadCursors() {
 // register our elements with skipper
 var elements = CR_getElements('*', 'cursor_stop');
 var len = elements.length;
 for (var i = 0; i < len; i++) {
 var element = elements[i]; 
 element.className = 'cursor_stop cursor_hidden';
 kibbles.skipper.append(element);
 }
 }
 function toggleComments() {
 CR_toggleCommentDisplay();
 reloadCursors();
 }
 function keysOnLoadHandler() {
 // setup skipper
 kibbles.skipper.addStopListener(
 kibbles.skipper.LISTENER_TYPE.PRE, updateCursor);
 // Set the 'offset' option to return the middle of the client area
 // an option can be a static value, or a callback
 kibbles.skipper.setOption('padding_top', 50);
 // Set the 'offset' option to return the middle of the client area
 // an option can be a static value, or a callback
 kibbles.skipper.setOption('padding_bottom', 100);
 // Register our keys
 kibbles.skipper.addFwdKey("n");
 kibbles.skipper.addRevKey("p");
 kibbles.keys.addKeyPressListener(
 'u', function() { window.location = detail_url; });
 kibbles.keys.addKeyPressListener(
 'r', function() { window.location = detail_url + '#publish'; });
 
 kibbles.keys.addKeyPressListener('j', gotoNextPage);
 kibbles.keys.addKeyPressListener('k', gotoPreviousPage);
 
 
 kibbles.keys.addKeyPressListener('h', toggleComments);
 
 }
 </script>
<script src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/code_review_scripts.js"></script>
<script type="text/javascript">
 function showPublishInstructions() {
 var element = document.getElementById('review_instr');
 if (element) {
 element.className = 'opened';
 }
 }
 var codereviews;
 function revsOnLoadHandler() {
 // register our source container with the commenting code
 var paths = {'svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b': '/scripts/DTI_analysis.py'}
 codereviews = CR_controller.setup(
 {"token":"3Cd3YLziNQwHJ6q0INBaXA2gZls:1366032649547","projectName":"pmx","domainName":null,"assetHostPath":"https://ssl.gstatic.com/codesite/ph","loggedInUserEmail":"vytautas.gapsys@gmail.com","profileUrl":"/u/110130407061490526737/","assetVersionPath":"https://ssl.gstatic.com/codesite/ph/14689258884487974863","projectHomeUrl":"/p/pmx","relativeBaseUrl":""}, '', 'svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b', paths,
 CR_BrowseIntegrationFactory);
 
 // register our source container with the commenting code
 // in this case we're registering the container and the revison
 // associated with the contianer which may be the primary revision
 // or may be a previous revision against which the primary revision
 // of the file is being compared.
 codereviews.registerSourceContainer(document.getElementById('lines'), 'svn82a17baf41be6e3dd11997ac8d7dff4272c3a37b');
 
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_DRAFT_PLATE, showPublishInstructions);
 
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_PUB_PLATE, pubRevealed);
 codereviews.registerActivityListener(CR_ActivityType.REVEAL_DRAFT_PLATE, draftRevealed);
 codereviews.registerActivityListener(CR_ActivityType.DISCARD_DRAFT_COMMENT, draftDestroyed);
 
 
 
 
 
 
 
 var initialized = true;
 reloadCursors();
 }
 window.onload = function() {keysOnLoadHandler(); revsOnLoadHandler();};

</script>
<script type="text/javascript" src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/dit_scripts.js"></script>

 
 
 
 <script type="text/javascript" src="https://ssl.gstatic.com/codesite/ph/14689258884487974863/js/ph_core.js"></script>
 
 
 
 
</div> 

<div id="footer" dir="ltr">
 <div class="text">
 <a href="/projecthosting/terms.html">Terms</a> -
 <a href="http://www.google.com/privacy.html">Privacy</a> -
 <a href="/p/support/">Project Hosting Help</a>
 </div>
</div>
 <div class="hostedBy" style="margin-top: -20px;">
 <span style="vertical-align: top;">Powered by <a href="http://code.google.com/projecthosting/">Google Project Hosting</a></span>
 </div>

 
 


 
 </body>
</html>

