<?php
$uploadpath="./uploads/";
if(isset($_POST['download'])){
	$filebasename=$_POST['filebasename'];
	if(strpos($_POST['download'],"Sequences")!==false){
		$outputfile=$filebasename."-Rotated.fasta";
		header("Content-type: text/x-fasta");
		header("Content-Disposition: attachment; filename=\"".$outputfile."\"");
		readfile($uploadpath.$outputfile);
		exit();
	}
	if(strpos($_POST['download'],"Blocks")!==false){
		$outputfile=$filebasename."-Blocks.csv";
		header("Content-type: text/csv");
		header("Content-Disposition: attachment; filename=\"".$outputfile."\"");
		readfile($uploadpath.$outputfile);
		exit();
	}
	if(strpos($_POST['download'],"Image")!==false){
		$outputfile=$filebasename."-Blocks.bmp";
		header("Content-type: image/bmp");
		header("Content-Disposition: attachment; filename=\"".$outputfile."\"");
		readfile($uploadpath.$outputfile);
		exit();
	}
}
if(isset($_POST['showlog'])){
	$outputfile="CSA_log.txt";
	header("Content-type: text/plain");
	header("Content-Disposition: inline; filename=\"".$outputfile."\"");
	readfile($uploadpath.$outputfile);
	exit();
}
header("Expires: Sat, 26 Jul 1997 05:00:00 GMT");
header("Last-Modified: ".gmdate("D, d M Y H:i:s")." GMT");
header("Cache-Control: no-store, no-cache, must-revalidate");
header("Cache-Control: pre-check=0, post-check=0, max-age=0");
header("Pragma: no-cache");
header("Expires: 0");
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Cyclic DNA Sequence Aligner</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link href="styles.css" rel="stylesheet" type="text/css" media="screen" />

<script type="text/javascript">
function toggleField(obj){
	var fileField=document.getElementById("file");
	var textField=document.getElementById("text");
	if(obj.value=="text"){
		fileField.style.display="none";
		textField.style.display="block";
		fileField.disabled=true;
		textField.disabled=false;
	}
	if(obj.value=="file"){
		fileField.style.display="block";
		textField.style.display="none";
		fileField.disabled=false;
		textField.disabled=true;
	}
}
function checkForm(){
	var field=null;
	var result=true;
	var minfield=document.getElementById("minblocksize");
	var maxfield=document.getElementById("maxblocksize");
	var minval=minfield.value;
	var maxval=maxfield.value;
	if(document.getElementById("typefile").checked)
		field=document.getElementById("file");
	if(document.getElementById("typetext").checked)
		field=document.getElementById("text");
	if(field.value.length==0){
		field.parentNode.className="error";
		result=false;
	} else
		field.parentNode.className="third";
	if(isNaN(minval) || minval.length==0 || minval<1)
		minfield.value=5;
	if(isNaN(maxval) || maxval.length==0 || (parseInt(maxval)<parseInt(minfield.value) && maxval!=0))
		maxfield.value=0;
	if(result==true){
		document.inputform.submitbutton.value="Please wait...";
		document.inputform.submitbutton.disabled=true;
		document.inputform.submit();
	}
	return result;
}
function goTo(startpos,endpos,seq){
	var table=document.getElementById("positionstable");
	var tableRows=table.tBodies[0].rows;
	var numRows=tableRows.length;
	var numCols=tableRows[0].cells.length;
	var cellHeight=(table.scrollHeight)/(numRows);
	var cellWidth=(table.scrollWidth)/(numCols-1);
	var val=0;
	var maxval=0;
	var minval=0;
	var indexminval=0;
	for(i=0;i<numRows;i++){ // get max row value
		if((val=parseInt(tableRows[i].cells[seq].innerHTML))>=maxval) maxval=val;
	}
	for(i=0;i<numRows;i++){ // slow sort rows
		minval=maxval;
		indexminval=0;
		for(j=0;j<(numRows-i);j++){
			if((val=parseInt(tableRows[j].cells[seq].innerHTML))<=minval){
				minval=val;
				indexminval=j;
			}
		}
		table.tBodies[0].appendChild(tableRows[indexminval]);
	}
	minval=endpos;
	indexminval=0;
	for(i=0;i<numRows;i++){ // get lowest row value higher than startpos
		if((val=parseInt(tableRows[i].cells[seq].innerHTML))>=startpos && val<=minval){
			minval=val;
			indexminval=i;
		}
	}
	document.getElementById("positionsframe").scrollTop=Math.round(indexminval*cellHeight);
	document.getElementById("positionsframe").scrollLeft=Math.round((seq-1)*cellWidth);
}
</script>

</head>

<?php
// If no uploaded file type is set, show the landing/submission page
if(!isset($_POST['type'])) {
	// print the request type and hostname to the log file
	$requestip=(getenv('HTTP_CLIENT_IP')?:getenv('HTTP_X_FORWARDED_FOR')?:getenv('HTTP_X_FORWARDED')?:getenv('HTTP_FORWARDED_FOR')?:getenv('HTTP_FORWARDED')?:getenv('REMOTE_ADDR'));
	$hostname=gethostbyaddr($requestip);
	$logfile=$uploadpath."CSA_log.txt";
	file_put_contents($logfile,date(DATE_RSS,$_SERVER['REQUEST_TIME'])."\t".$hostname."\n",FILE_APPEND);
?>

<body onload="document.getElementById('text').style.display='none';">

<div id="content">
<div id="back">

<div id="header">
	<div id="menu">
		<ul>
			<li id="button1"><a href="#description">Home</a></li>
			<li id="button2"><a href="#submition">Submit</a></li>
			<li id="button3"><a href="#document">Documentation</a></li>
			<li id="button4"><a href="#download">Downloads</a></li>
			<li id="button5"><a href="#reference">Reference</a></li>
			<li id="button6"><a href="#contact">Contact us</a></li>
		</ul>
	</div>
	<a href="http://www.ipatimup.pt/"><img src="ipatimup-logo.gif" style="border:0px;height:5%;width:15%;position:absolute;top:10%;right:22%;z-index:5;"/></a>
	<a href="http://kdbio.inesc-id.pt/"><img src="kdbio-logo.gif" style="border:0px;height:10%;width:15%;position:absolute;top:8%;right:42%;z-index:5;"/></a>
	<div id="logo">
		<form id="showlog" method="post" action="<?php print($_SERVER['PHP_SELF']); ?>"><input type="hidden" name="showlog" value="1"/></form>
		<h1>
		<div onclick="javascript:document.getElementById('showlog').submit();" style="text-shadow: 0px 0px 10px black;">
		CSA
		</div>
		</h1>
		<h2>Cyclic DNA Sequence Aligner</h2>
	</div>
</div>

<div id="main">
 <div id="right">
	<h3>Multiple Sequence Alignment Tools</h3>
	<ul><li>
		<h4><a target="_blank" href="http://www.clustal.org/clustal2/">ClustalW</a></h4><br/>
		<h4><a target="_blank" href="http://darlinglab.org/mauve/">MAUVE</a></h4><br/>
		<h4><a target="_blank" href="http://alggen.lsi.upc.es/recerca/align/mgcat/">M-GCAT</a></h4><br/>
		<h4><a target="_blank" href="http://bio.math.berkeley.edu/mavid/download/">MAVID</a></h4><br/>
		<h4><a target="_blank" href="http://toolcoconut.org/">CoCoNUT</a></h4><br/>
		<h4><a target="_blank" href="http://biocomp.iis.sinica.edu.tw/new/sinicview.php">SinicView</a></h4><br/>
	</li></ul>
	<h3>Sequence Genome Viewer</h3>
	<ul><li>
		<h4><a target="_blank" href="http://www.sanger.ac.uk/Software/Artemis/circular/">DNAPlotter</a></h4><br/>
	</li></ul>
 </div>
 <div id="left">
	<h2>Cyclic DNA Sequence Aligner</h2><br/>
	
	<h4><a name="description">DESCRIPTION</a></h4><br/>
	This tool finds the optimal rotation for a set of <b>circular DNA sequences</b> that are going to be aligned.
	It is very well suited to apply to <i>mitochondrial genome alignments</i>,
	and also for plasmids, chloroplasts, circular viruses and bacterial chromosomes.
	<br/><br/>
	The best rotation is calculated based on the longest chain of non-repeated blocks that belongs to all the sequences.
	These maximum common blocks are obtained with the help of a <i>generalized cyclic suffix tree</i> data structure.
	<br/><br/>
	As of yet, the current version of this tool does not perform the multiple sequence alignment itself.
	But some suggestions of other highly efficient <i>multiple sequence alignment tools</i> to perform this task are presented on the right.
	<br/><br/><br/>

	<h4><a name="submition">SUBMISSION</a></h4><br/>
	The sequences can be submitted in Multi-FASTA format by uploading a file or by copying/pasting the sequences in the text window.
	There is a limit of 64 sequences and 5 MB per file.
	<!-- The length of the blocks displayed in the output can be specified by the user. -->
	<br/><br/>

	<form name="inputform" enctype="multipart/form-data" method="post" action="<?php print($_SERVER['PHP_SELF']); ?>">
	<input type="hidden" name="MAX_FILE_SIZE" value="5000000" />
	
	<table class="inputtable" id="formtable">
	<tr>
	<td class="first">Sequences</td>
	<td class="third" colspan="2">
	<input type="file" id="file" name="file" size="45" />
	<textarea id="text" name="text" wrap="virtual" rows="5" cols="45" style="word-wrap:break-word;"></textarea>
	</td>
	</tr>
	<tr>
	<td class="first">Options</td>
	<td class="second">
	<input type="radio" value="file" name="type" id="typefile" checked="checked" onclick="toggleField(this);" />
	<label for="typefile">Upload File</label>
	<br/>
	<input type="radio" value="text" name="type" id="typetext" onclick="toggleField(this);" />
	<label for="typetext">Upload Text</label>
	</td>
	<td class="third">
	<input type="text" id="minblocksize" name="minblocksize" size="2" maxlength="4" value="10" readonly="readonly" />
	<label for="minblocksize">Minimum block size (<i>&ge;1</i>)</label>
	<br/>
	<input type="text" id="maxblocksize" name="maxblocksize" size="2" maxlength="4" value="0" readonly="readonly" />
	<label for="maxblocksize">Maximum block size (<i>0=&infin;</i>)</label>
	</td>
	</tr>
	<tr>
	<td class="third" colspan="3">
	<input type="button" name="submitbutton" value="Submit" onclick="return checkForm();" />
	</td>
	</tr>
	</table>

	</form>

	<br/><br/>

	<h4><a name="document">DOCUMENTS</a></h4><br/>
	<a href="UserManual.pdf">User Manual</a>
<!--	<br/><a href="SupplementaryMaterial.zip">Supplementary Material</a> -->
	<br/><br/><br/>

	<h4><a name="download">DOWNLOADS</a></h4><br/>
	<a href="SourceCode.zip">Source code</a>
	<br/><a href="Examples.zip">Example sequences</a>
	<br/>Use '<i>unzip</i>' to unpack and '<i>make</i>' to compile.
	<br/>Run with: <i>CSA R sequences.fasta</i>
	<br/><br/><br/>

	<h4><a name="reference">REFERENCE</a></h4><br/>
	Fernandes F., Pereira L., Freitas A.T.:
	<br/><a href="http://www.biomedcentral.com/1471-2105/10/230">CSA - An efficient algorithm to improve circular DNA multiple alignment.</a>
	<br/>BMC Bioinformatics, 2009, 10:230.
	<br/><br/><br/>

	<h4><a name="contact">CONTACT</a></h4><br/>
	Please send your comments, suggestions, bug reports or questions to:
	<br/><a href="mailto:csatool@kdbio.inesc-id.pt">csatool@kdbio.inesc-id.pt</a>
	<br/><br/><br/>

 </div>
</div>

<div id="footer">
<p>Last updated: <?php print(date("d M Y",filemtime(__FILE__))); ?></p>
<!--
<p>Copyright &copy; 2008. <a href="#">Privacy Policy</a> | <a href="#">Terms of Use</a> | <a href="http://validator.w3.org/check/referer" title="This page validates as XHTML 1.0 Transitional"><abbr title="eXtensible HyperText Markup Language">XHTML</abbr></a> | <a href="http://jigsaw.w3.org/css-validator/check/referer" title="This page validates as CSS"><abbr title="Cascading Style Sheets">CSS</abbr></a></p>
	<p>Design by <a href="http://www.metamorphozis.com/" title="Flash Templates">Flash Templates</a>, coded by <a href="http://www.flashtemplatesdesign.com" title="Free Flash Templates">Free Flash Templates</a>
		</p>
-->
</div>
<br/>
</div>
</div>
<br/>
<!--
<div style="text-align: center; font-size: 0.75em;">Design downloaded from <a href="http://www.freewebtemplates.com/">free website templates</a>.</div></body>
-->
</body>
</html>

<?php
// Else, if an uploaded file type is set, run the program and show the results page
} else {
	// Clean up old uploaded/created files
	//echo "<div id=\"codetext\"><pre>\n";
	$fileslist=scandir($uploadpath);
	foreach($fileslist as $file){
		if($file=="." || $file=="..") continue;
		$file=$uploadpath.$file;
		// delete files larger than 5MB
		if( floor(filesize($file)/(1024*1024)) > 5 ){
			unlink($file);
			continue;
		}
		if($file==$uploadpath."CSA_log.txt") continue;
		// delete files older than 48 hours
		if( floor((time()-filemtime($file))/(60*60)) > 48 ){
			unlink($file);
			continue;
		}
		//echo floor((time()-filemtime($file))/(60*60)) . "h\t";
		//echo floor(filesize($file)/(1024*1)) . "KB\t";
		//echo "$file\t".date(filemtime($file))."\n";
	}
	//echo "</pre></div>\n";
?>

<body>
<div id="content">
<div id="back">
<div id="main">
<?php
// file_exists() ; is_writable()
if($_POST['type']=="file"){
	$filename=basename($_FILES['file']['name']);
	$filename=preg_replace("/[^a-zA-Z0-9-._]/","",$filename);
	$filesource=$_FILES['file']['tmp_name'];
	$filepath=$uploadpath.$filename;
	move_uploaded_file($filesource,$filepath);
}
if($_POST['type']=="text"){
	$filepath=$uploadpath."sequence".time().".txt";
	$file=fopen($filepath,"w");
	fwrite($file,$_POST['text']);
	fclose($file);
}
$fileext=pathinfo($filepath,PATHINFO_EXTENSION);
$filebasename=basename($filepath,".".$fileext);
//$filepath=$uploadpath.$filename;
//$commandargs=$filepath." -S ".$_POST['minblocksize']." -W ".$_POST['maxblocksize'];
$commandargs=$filepath;
$outputimage=$uploadpath.$filebasename."-Blocks.bmp";
$mapfilename=$uploadpath.$filebasename."-imagemap.txt";
$datafilename=$uploadpath.$filebasename."-positions.txt";
if(file_exists($outputimage)) unlink($outputimage);
if(file_exists($mapfilename)) unlink($mapfilename);
if(file_exists($datafilename)) unlink($datafilename);
?>
<div id="codetext"><pre>
<?php
$cmdexec="timeout -s 9 1h ./CSA R";
if(stripos($_SERVER["SERVER_SOFTWARE"],"Win")) $cmdexec="Debug\CyclicSequenceAligner.exe";
$cmdexec.=" ".$commandargs;
$command=popen($cmdexec,"r");
if($command!==false) {
	ob_start();
	$skipchar=0;
	while( ($char=fgetc($command)) !== false ) {
		if(ord($char)==27) $skipchar=10; // ESC char
		$char=htmlspecialchars($char);
		if($char=="\n"){
			$char="\r\n";
			$skipchar=0;
		}
		if($skipchar<=0) print($char);
		$skipchar--;
		ob_flush();
		flush();
		usleep(1000);
	}
	ob_end_flush();
	pclose($command);
} else print("Error while executing command \"".$cmdexec."\"");
unlink($filepath);
?>
</pre></div>
<br/>
<table id="imagestable" class="resultstable">
<tr><td>
<?php
//$mapfilename="imagemap.map";
$mapfile=fopen($mapfilename,"r");
//$mapfile=fopen($mapfilename,"w");
//if($datafile && $mapfile){
if($mapfile){
	$nx=0;
	$ny=0;
	fscanf($mapfile,"%d %d\n",$nx,$ny);
	$ximg=explode(" ",trim(fgets($mapfile)));
	$xseq=explode(" ",trim(fgets($mapfile)));
	$yup=explode(" ",trim(fgets($mapfile)));
	$ydown=explode(" ",trim(fgets($mapfile)));
	printf("<map name=\"imagemap\" id=\"imagemap\">\n");
	for($j=0;$j<$ny;$j++){
		for($i=0;$i<($nx-1);$i++){
			$coords=($ximg[$i]).",".($yup[$j]).",".($ximg[$i+1]-1).",".($ydown[$j]);
			$href="javascript:goTo(".($xseq[$i]).",".($xseq[$i+1]).",".($j+1).");";
			printf("<area shape=\"rect\" coords=\"".$coords."\" href=\"".$href."\"/>\n");
		}
	}
	printf("</map>\n");
	fclose($mapfile);
}
list($imagewidth,$imageheight)=getimagesize($outputimage);
$frameheight=($imageheight+4);
?>
<img class="imageframe" src="<?php print($outputimage."?".time()); ?>" usemap="#imagemap" alt="Aligned blocks" />
</td>
<td>
<div id="positionsframe" class="positionsframe">
<table id="positionstable" class="positionstable">
<?php
$datafile=fopen($datafilename,"r");
if($datafile){
	$ns=0;
	fscanf($datafile,"%d\n",$ns);
	print("<thead><tr><th>Size</th>");
	for($i=1;$i<=$ns;$i++) print("<th>Seq.&nbsp;".$i."</th>");
	print("</tr></thead>\n<tbody>\n");
	while(!feof($datafile)){
		$dataline=explode(" ",trim(fgets($datafile)));
		if(count($dataline)<4) break;
		$colorval="#";
		$colorval.=str_pad(dechex($dataline[0]),2,"0",STR_PAD_LEFT);
		$colorval.=str_pad(dechex($dataline[1]),2,"0",STR_PAD_LEFT);
		$colorval.=str_pad(dechex($dataline[2]),2,"0",STR_PAD_LEFT);
		$coloradd=255-max($dataline[0],$dataline[1],$dataline[2])/2;
		$lowcolorval="#";
		$lowcolorval.=str_pad(dechex(($dataline[0]/2)+$coloradd),2,"0",STR_PAD_LEFT);
		$lowcolorval.=str_pad(dechex(($dataline[1]/2)+$coloradd),2,"0",STR_PAD_LEFT);
		$lowcolorval.=str_pad(dechex(($dataline[2]/2)+$coloradd),2,"0",STR_PAD_LEFT);
		print("<tr bgcolor=\"".$lowcolorval."\"><td bgcolor=\"".$colorval."\">".$dataline[3]."</td>");
		for($i=1;$i<=$ns;$i++) print("<td>".$dataline[(3+$i)]."</td>");
		print("</tr>\n");
	}
	print("</tbody>\n");
	fclose($datafile);
}
?>
</table>
</div>
</td></tr>
</table>
<table id="buttonstable" class="resultstable"><tr>
<form name="resultsform" method="post" action="<?php print($_SERVER['PHP_SELF']); ?>">
	<input type="hidden" name="filebasename" value="<?php print($filebasename); ?>" />
<td><input type="submit" name="download" value="Download Rotated Sequences" id="downloadbutton1" /></td>
<td><input type="submit" name="download" value="Download Original Blocks Info" id="downloadbutton2" /></td>
<td><input type="submit" name="download" value="Download Image" id="downloadbutton3" /></td>
<td><input type="submit" value="Go Back" /></td>
</form>
</tr></table>
<br/>
<script type="text/javascript">
var outputtext=document.getElementById("codetext").innerHTML;
if(outputtext.indexOf("Done!")==-1 || outputtext.indexOf("ERROR")!=-1 || outputtext.indexOf("USAGE")!=-1) {
	document.getElementById("downloadbutton1").disabled=true;
	document.getElementById("downloadbutton2").disabled=true;
	document.getElementById("downloadbutton3").disabled=true;
	document.getElementById("imagestable").className="hidden";
	document.getElementById("positionstable").className="hidden";
	document.getElementById("linkstable").className="hidden";
}
document.getElementById("positionsframe").style.height=<?php print($frameheight); ?>;
<!--
document.getElementById("darkbackground").style.width=document.body.scrollWidth;
document.getElementById("darkbackground").style.height=document.body.scrollHeight;
document.getElementById("waitmessage").style.top=Math.round(document.body.scrollTop + document.body.clientHeight/2);
document.getElementById("waitmessage").style.left=Math.round(document.body.scrollLeft + document.body.clientWidth/2);
-->
</script>

</div>
<div id="footer"></div>
<br/>
</div>
</div>
<br/>
</body>
</html>

<?php } ?>
