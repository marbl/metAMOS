import sys, os, string
ROOT = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, os.path.join(ROOT, '..'))
#sys.path.append(ROOT+"/lib")
import  markup, datetime
from pygooglechart import StackedVerticalBarChart
from pygooglechart import PieChart2D
from pygooglechart import PieChart3D
from pygooglechart import Axis
from pygooglechart import StackedHorizontalBarChart, StackedVerticalBarChart, \
         GroupedHorizontalBarChart, GroupedVerticalBarChart
from datetime import datetime, date, time


import settings
import helper

def outputLibraryInfo(html_prefix, markupFile, outputHeader, libcnt, format, mated, interleaved, mmin, mmax, outputFastQC):
   if outputHeader:
      markupFile.table(border="1")
      markupFile.tr()
      markupFile.th("Library #")
      markupFile.th("Format")
      markupFile.th("Mated?")
      markupFile.th("Min Insert Size")
      markupFile.th("Max Insert Size")
      if outputFastQC:
         markupFile.th("First FastQC Report")
         markupFile.th("Second FastQC Report") 
      markupFile.tr.close()

   markupFile.tr()
   markupFile.add("<td align=\"left\">%d</td>"%(libcnt))
   markupFile.add("<td align=\"left\">%s</td>"%(format))
   markupFile.add("<td align=\"left\">%s</td>"%(mated))
   markupFile.add("<td align=\"right\">%d</td>"%(mmin))
   markupFile.add("<td align=\"right\">%d</td>"%(mmax))

   if outputFastQC and os.path.exists("lib%d.1.fastqc/fastqc_report.html"%(libcnt)):
      if mated.lower() == "true":
         if interleaved.lower() == "true":
            markupFile.td('<a target="_blank" href="lib%d.1.fastqc/fastqc_report.html">interleaved</a>'%(libcnt))
            markupFile.td("NA")
         else:
            markupFile.td('<a target="_blank" href="lib%d.1.fastqc/fastqc_report.html">left</a>'%(libcnt))
            markupFile.td('<a target="_blank" href="lib%d.2.fastqc/fastqc_report.html">right</a>'%(libcnt))
      else:
         markupFile.td('<a target="_blank" href="lib%d.1.fastqc/fastqc_report.html">unmated</a>'%(libcnt))
         markupFile.td("NA")
   markupFile.tr.close()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "usage: create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm> <Utils dir> <run dir> <# of libs> <taxa level of classifications>"
        sys.exit(0)
    rund = sys.argv[6]
    utils = sys.argv[5]

    prefix = sys.argv[3]
    html_prefix = prefix
    prefix = prefix.replace("/html/", "")
    MA_dir = prefix
    MA_dir = MA_dir.replace("/Postprocess/out","")
    ref_asm = sys.argv[4]
    mp = open(sys.argv[1],'r')
    nLibs = int(sys.argv[7])
    taxa_level = sys.argv[8]
    #mp2 = open(sys.argv[1].replace("s12","s3"),'r')    

    # set working dir
    os.chdir(html_prefix)

    if not os.path.exists(html_prefix+"asmstats.out"):
        libPath = rund.replace("bin", "lib")
        os.system("perl -I %s %s/perl/statistics.pl %s > %sasmstats.out"%(libPath,utils,sys.argv[4],html_prefix))        
    report = open(html_prefix+"asmstats.out",'r')

    # define colors and style elements
    textColor = "solid black"
    backgroundColor = "#E8E8E8"
    selectedColor = "yellow"
    mouseOverColor = "#cdcdcd"

    initialStep = "Annotate"

    steps = []
    steps.append("Preprocess")
    steps.append("Assemble")
    steps.append("MapReads")
    steps.append("FindORFS")
    steps.append("FindRepeats")
    steps.append("Scaffold")
    steps.append("FindScaffoldORFS")
    steps.append("Abundance")
    steps.append("Annotate")
    steps.append("Propagate")
    steps.append("Classify")

    step_status = {}
    step_status["Preprocess"] = "OK"
    step_status["Assemble"] = "OK"
    step_status["MapReads"] = "OK"
    step_status["FindORFS"] = "OK"
    step_status["FindRepeats"] = "OK"
    step_status["Scaffold"] = "OK"
    step_status["FindScaffoldORFS"] = "OK"
    step_status["Abundance"] = "FAIL"
    step_status["Annotate"] = "FAIL" 
    step_status["Propagate"] = "NONE"
    step_status["Classify"] = "NONE"
    ##get status of each step from Log dir
    for step in steps:
        stepn = step.lower()
        started = False
        completed = False
        skipped = False
        if os.path.exists("%s/Logs/%s.started"%(MA_dir,stepn)):
            started = True
        if os.path.exists("%s/Logs/%s.ok"%(MA_dir,stepn)):
            completed = True
        if os.path.exists("%s/Logs/%s.skip"%(MA_dir,stepn)):
            skipped = True
        if started and completed and not skipped:
            step_status[step] = "OK"
        elif started and completed and skipped:
            step_status[step] = "SKIP"
        elif started and not completed:
            step_status[step] = "FAIL"
        elif not started:
            step_status[step] = "NONE"
        else:
            step_status[step] = "NONE" 

    #+male1  /cbcb/project-scratch/sergek/metAMOS/individualAsms/m1_asm      proba   b-      metaphyler=1
    ## call Dan's script, for now on a single sample/run
    cpfile = open("%s/plot.tab"%(html_prefix),'w')
    cpfile.write("+sample1\t%s\tproba\tb-\tmetaphyler=1\n"%(MA_dir))
    cpfile.close()
    os.system("python %s/python/create_plots.py %s/plot.tab proba1"%(utils,html_prefix))

    ##update counts
    #count reads
    os.system("grep -c \">\" %s/Preprocess/out/*.fasta > readcount.txt"%(MA_dir))
    readcount = open("readcount.txt",'r').read().replace("\n","")  
    #print readcount
    os.system("rm readcount.txt")
    #count contigs
    os.system("grep -c \">\" %s/Assemble/out/proba.asm.contig > contigcount.txt"%(MA_dir))
    contigcount = open("contigcount.txt",'r').read().replace("\n","")  
    #print contigcount
    os.system("rm contigcount.txt")
    #count scaffolds
    os.system("grep -c \">\" %s/Scaffold/out/proba.scaffolds.final > scafcount.txt"%(MA_dir))
    scaffoldcount = open("scafcount.txt",'r').read().replace("\n","")  
    #print scaffoldcount
    os.system("rm scafcount.txt")
    #count scaffolds
    os.system("grep -c \">\" %s/Scaffold/out/proba.motifs > motifcount.txt"%(MA_dir))
    motifcount = open("motifcount.txt",'r').read().replace("\n","")  
    #print motifcount
    os.system("rm motifcount.txt")
    #count ORFs
    os.system("grep -c \">\" %s/FindORFS/out/proba.orfs.faa > orfcount.txt"%(MA_dir))
    orfcount = open("orfcount.txt",'r').read().replace("\n","")  
    #print orfcount
    os.system("rm orfcount.txt")
    ##copy stuff
    for step in steps:
        step = step.lower()
        os.system("cp %s/javascript/%s.js %s/."%(utils,step,html_prefix))
    os.system("cp %s/Logs/COMMANDS.log %s/pipeline.commands"%(MA_dir,html_prefix))
    os.system("cp %s/pipeline.run %s/pipeline.summary"%(MA_dir,html_prefix))

    # generate dynamic java scripts
    # first classify and propagate
    os.system("python %s/python/get_classify_stats.py %s/propagate.in.clusters %s/propagate.out.clusters %s/DB/tax_key.tab %s classify.js propagate.js %s"%(utils, html_prefix, html_prefix, utils, html_prefix, taxa_level)) 

    # generate preprocess
    preprocess = markup.page()
    preprocess.p()

    nQC = 0
    for i in range(1, nLibs + 1):
        if os.path.exists("%s/lib%d.1.fastqc/fastqc_report.html"%(html_prefix, i)):
            nQC = nQC + 1

    summary = open("%s/pipeline.ini"%(MA_dir), 'r')
    libcnt = 1
    format = None 
    mmin = 0
    mmax = 0
    mated = False
    interleaved = False
    innie = False
    linkerType = ""
    libadded = False
    firstLib = True
    for line in summary:
       line = line.replace("\n","")
       if "#" in line:
          continue
       elif "asmcontigs:" in line:
          asmc = line.replace("\n","").split("\t")[-1]
          if len(asmc) <= 2:
             continue
          preprocess.add("Supplied assembly: %s"%(asmc))
          preprocess.br()
       elif "format:" in line:
          if format and not libadded:
             outputLibraryInfo(html_prefix, preprocess, firstLib, libcnt, format, mated, interleaved, mmin, mmax, nQC > 0)
             libcnt += 1
          libadded = False
          format = line.replace("\n","").split("\t")[-1]
       elif "mated:" in line:
          mated = line.replace("\n","").split("\t")[-1]
       elif "innie:" in line:
         innie = line.replace("\n","").split("\t")[-1]
       elif "linker:" in line:
         linkerType = line.replace("\n","").split("\t")[-1]
       elif "interleaved:" in line:
         interleaved = line.replace("\n","").split("\t")[-1]
       elif "f1:" in line:
          data = line.split("\t")

          inf = data[1].split(",")
          mmin = int(inf[1])
          mmax = int(inf[2])
       elif "f2:" in line:
          data = line.split("\t")

          inf = data[1].split(",")
          mmin = int(inf[1])
          mmax = int(inf[2])

          libadded = True
       elif "frg" in line:
          data = line.split("\t")
          mated = False
          libadded = True
    if format and not libadded:
       outputLibraryInfo(html_prefix, preprocess, firstLib, libcnt, format, mated, interleaved, mmin, mmax, nQC > 0)
    preprocess.table.close()
    summary.close()

    preprocess_out = open("%s/preprocess.js"%(html_prefix), 'w')
    preprocess_out.write("preprocessHTML = '%s'"%(preprocess.__str__().replace("\n", "\\\n")))
    preprocess_out.close()
    
    # todo, need to add report for MapReads including # reads mapped (%), contig coverage histogram, and % reads between contigs and number of links histogram. Also re-estimated insert sizes for each lib
    mapreads = markup.page()
    mapreads.img(src_="hist_ctgcvg.png",height_="100%",width_="100%")
    mapreads_out = open("%s/mapreads.js"%(html_prefix), 'w')
    mapreads_out.write("mapreadsHTML = '%s'"%(mapreads.__str__().replace("\n", "\\\n")))

    ##This will create ScaffoldSizes.png,ContigSizes.png
    
    ##create code to automatically generate .js files for HTML report
    ## let's start with Abundance
    #abundance_js = open(html_prefix+"abundance.js",'w')
    
    rdata = []
    for line in report:
        rdata.append(line)
       
    if not os.path.exists(html_prefix+"covstats.out"):
        os.system("%s/analyze-read-depth -x 2 %s > %scovstats.out"%(rund,sys.argv[2],html_prefix))
    ff = open(html_prefix+"covstats.out",'r')
    covdata = []
    #covdata = ff.readlines()
    #zflag = 0
    for line in ff:
        covdata.append(line)
    
    if not os.path.exists(html_prefix+"stats.out"):
        os.system("%s/astats %s > %sstats.out"%(rund,sys.argv[2],html_prefix))
    dd = open(html_prefix+"stats.out",'r')
    ddata = dd.readlines()

    ddf = ""
    bflag = 0
    cflag = 0
    rflag = 0
    for line in ddata:
         if "<body>" in line:
             bflag = 1
             ddf += "<table border=\"1\">"
         elif "[Contigs]" in line:
             cflag = 1
             ddf += line
         elif "[Small" in line:
             cflag = 0
         elif "[Big" in line:
             cflag = 0
         elif "[Reads]" in line:
             rflag = 1
             ddf += line
         elif cflag != 0 or rflag != 0:
             if "N50" not in line:
                 ddf += line
         else:
             continue
    mp.readline()
    cov = []
    abund = []
    ids = []
    mpd = mp.readlines()
    phylum = False
    for line in mpd:
        if ">phylum" in line:
            phylum = True
        if len(line) < 3:
            continue
        data = line.split("\t")
        if len(data) < 2:
            continue
        if not phylum:
            continue
        ids.append(data[0])
        #cov.append(int(float(data[2])))
        abund.append(100*float(data[1]))


        
    # Create a chart object of 200x100 pixels
#    chart2 = StackedVerticalBarChart(600, 300)
    chart2 = PieChart3D(550, 300)
    
    chart2.set_colours(["0A8C8A","EBB671","DE091A"])
    # Add some data
    chart2.add_data(abund)
#    chart2.add_data(abund2)    

    # Assign the labels to the pie data
    chart2.set_pie_labels(ids)
#    chart2.set_axis_labels(Axis.BOTTOM, ids)

    # Download the chart
    try:
        chart2.download(html_prefix+'abund.png')
    except:
        print "Warning: could not download abund.png"

    chart = GroupedHorizontalBarChart(600, 500, x_range=(0,100))
    chart.set_bar_width(30)
    chart.set_colours(("FF0000","0000FF"))
    #chart.set_colours_within_series(("FF0000","FF7F00","FFFF00","00FF00","0000FF","8B00FF","FFFFFF"))

#    chart.set_colours(['00ff00', 'ff0000'])
    chart.add_data(abund)
    #chart.add_data(abund2)
    #chart.add_data([1,4,9,16,25])
    #category = ["even","staggered"]
    #chart.set_legend(category)            
    ids = ids[::-1]

    chart.set_axis_labels(Axis.LEFT, ids)
    chart.set_axis_labels(Axis.BOTTOM, [0,10,20,30,40,50,60,70,80,90,100])    
    #chart.set_axis_labels(Axis.LEFT, ["0)
    index = chart.set_axis_labels(Axis.BOTTOM, ['Abundance (%)'])
    chart.set_axis_style(index, '202020', font_size=10, alignment=0)
    chart.set_axis_positions(index, [50])

    index = chart.set_axis_labels(Axis.LEFT, ['Phylum'])
    chart.set_axis_style(index, '202020', font_size=10, alignment=0)
    chart.set_axis_positions(index, [50])
    
    try:
        chart.download(html_prefix+'bar-phylum.png')
    except:
        print "Warning: could not download bar-phylum.png"
    
    dt = datetime.now()
    ds = dt.strftime("%A, %d. %B %Y %I:%M%p")
    title = "metAMOS: a metagenomic assembly pipeline for AMOS"
    # header is blank for now
    header = []

    # write the javascript we need on the page
    script = []
    script.append("<script type=\"text/javascript\">")
    script.append("moverColor = \"%s\""%(mouseOverColor))
    script.append("selectedColor = \"%s\""%(selectedColor))
    script.append("defaultColor = \"%s\""%(backgroundColor))
    script.append("bgcolor = \"white\"")
    script.append("selected = \"\"")
    script.append("")
    script.append("function load(fileName) {")
    script.append("   if (selected != \"\" && selected == fileName) {")
    script.append("      return;");
    script.append("   }");
    script.append("   $('#krona').html(window[fileName.concat('HTML')]);")
    script.append("   $('#links tr').each(function(){")
    script.append("      $(this).find('td').each(function(){")
    script.append("      if (($(this).attr('id')).toLowerCase() == fileName) {")
    script.append("         $(this).css(\"background-color\", selectedColor);")
    script.append("         bgcolor = selectedColor;")
    script.append("         selected = fileName;")
    script.append("      } else {")
    script.append("         $(this).css(\"background-color\", defaultColor);")
    script.append("      }")
    script.append("   })")
    script.append("})")
    script.append("}")
    script.append("")
    script.append("function mover(aa) {")
    script.append("   bgcolor = aa.style.backgroundColor;")
    script.append("   aa.style.backgroundColor = moverColor;")
    script.append("}")
    script.append("")
    script.append("function mout(aa) {")
    script.append("   if (bgcolor == selectedColor) {")
    script.append("      aa.style.backgroundColor = selectedColor;")
    script.append("   } else {")
    script.append("      aa.style.backgroundColor = defaultColor;")
    script.append("   }")
    script.append("}")
    script.append("</script>")

    #<link rel="shortcut icon" href="../assets/ico/favicon.ico">
    #<link rel="apple-touch-icon-precomposed" sizes="144x144" href="../assets/ico/apple-touch-icon-144-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" sizes="114x114" href="../assets/ico/apple-touch-icon-114-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" sizes="72x72" href="../assets/ico/apple-touch-icon-72-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" href="../assets/ico/apple-touch-icon-57-precomposed.png">

    # generate the dictionary of javascript pages we need
    scripts = {}
    scripts["file://%s/javascript/jquery-latest.js"%(utils)] = "javascript"
    scripts["http://code.jquery.com/jquery-latest.js"] = "javascript"
    # now a javascript for each page
    for step in steps:
       scripts["%s.js"%(step.lower())] = "javascript"

    # attributes for the body tax
    body = {}
    body["onload"] = "load('%s')"%(initialStep.lower())

    # the footer for the page

    footer = ""#"Generated %s"%(ds)
    styles = ( 'style2.css')#'./html/bootstrap.css', './html/boostrap-responsive.css')#'style2.css')#'layout.css', 'alt.css', 'images.css' )
    #styles = ( 'layout.css', 'alt.css', 'images.css' )
    #meta = ('viewport':"width=device-width, initial-scale=1.0",'description':'','author':'')

    # create HTML page
    page = markup.page( )
    # add the javascript free form
    page.add("\n".join(script))
    # initialize the rest of the body/html headers
    page.init( title=title,   \
               script=scripts,\
               header='\n'.join(header), \
               bodyattrs=body, \
               footer=footer )
    
#    page.br()
    #page.div( class_ = 'navbar navbar-fixed-top')
    #page.div( class_ = 'navbar-inner')
    page.div( class_ = 'container', style_='width:100%;height=100%')
    #page.div.close()
    #page.div.close()
    #page.div.close()
    #page.div( id_='page' )
    page.div( id_='header', style="background-color:#B8B8B8;clear:both;text-align:center;width:100%;height:7%;border:1px solid black") 
    page.h1("<u>MetAMOS <font color=\"blue\">v1.0</font> metagenomic assembly & analysis report</u>" , style_="text-align:center;vertical-align:top")
    page.h1.close()
    #page.div( id_='title')# style="vertical-align:bottom;")
    #page.font( size=14)
    #page.br( )
    page.div.close()

    #<frameset rows="40%,60%"cols="80%,20%">
    #<frame src="report.krona.html">
    #<frameset rows="20%,20%">
    #<frame src="report.krona.html">
    #<frame src="report.krona.html">
    #</frameset>
    #<frameset rows="60%" cols="40%,40%">
    #<frame src="report.krona.html">
    #<frame src="report.krona2.html">
    #</frameset>
    #<frameset rows="20%,20%,20%" cols="20%">
    #<frame src="report.krona.html">
    #<frame src="report.krona.html">
    #<frame src="report.krona.html">
    #</frameset>
    #</frameset>
    #</html>    
    #font, estimate
    #page.div.close()
    #page.frameset( rows_="40%,60%", cols_ = "80%,20%")
    #page.div()
    page.div( id_="menu", style_="background-color:%s;text-align:left;width:10%%;height:88%%;float:left;border:1px %s"%(backgroundColor, textColor))
    
    #items = ["<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>", ]
    #page.ul()
    page.p("<u>Pipeline status</u><br>")
    page.table( id_="links" )
    for step in steps:
       page.tr()
       status = "NA"
       color = "gray"
       try:
          if step_status[step] == "OK":
             status = "OK"
             color = "green"
          elif step_status[step] == "FAIL":
             status = "FAIL"
             color = "red"
          elif step_status[step] == "SKIP":
             status = "SKIPPED"
             color = "blue"
       except KeyError:
          continue

       tableHTML = []
       tableHTML.append("<td id=\"%s\" onmouseover=\"mover(this);\" onmouseout=\"mout(this);\">"%(step.lower()))
       tableHTML.append("<a href=\"javascript:void(0)\" onclick=\"load('%s');\">%s</a>"%(step.lower(), step))
       tableHTML.append("<br>")
       tableHTML.append("<font color=\"%s\">%s</font>"%(color, status))
       tableHTML.append("</td>")
       page.add("\n".join(tableHTML))
       page.tr.close()

    page.table.close()
    #page.ul.close()
    page.div.close()
    #page.li("<a target=\"_blank\" href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>")    
    #page.li("<a href=\"http://cbcb.umd.edu/~mpop/Software.shtml\">Related software</a>")
    #page.li("<a href=\"http://cbcb.umd.edu\">CBCB</a>")
    
#    page.li( items[1] )
    #page.ul.close( )
    #page.div.close()
    page.div( id_="content", style="background-color:#FFFFFF;float:left;width:58%%;height:12%%;border:1px %s"%(textColor))

    # TODO: do we want this? also, test -BDO
    #if os.path.exists("%s/Annotate/out/report.krona.html"%prefix):
    #    page.iframe(src="%s/Annotate/out/report.krona.html"%prefix, width="100%", height="600px")

    #page.div( id_='wrapper')
    #page.div( id_="content")
    
    #page.a( "Reference assembly:", class_='internal', href='%s'%(ref_asm) )          

    #paragraphs = ( ddf )
    cnt = 0
    #table_html = ""
    if 0:
        page.table(border="1",style_="width=80%")
        for contig in covdata:
            if cnt == 0:
                #table_html += "<table border=\"1\">\n"
                page.tr()
                page.td("High Coverage Contig ID")
                page.td("Coverage")
                page.tr.close()
                #table_html += "<tr><td> High Coverage Contig ID </td> <td> Coverage </td></tr>\n"
            page.tr()
            page.td(contig.split("\t")[0])
            page.td(contig.split(" ")[-1])


            cnt +=1

    #page.table.close()        
    #page.p( paragraphs )
    #page.p(style_="font-size:6px")
    page.table(border="2",width="80%")
    #page.tr()
    #page.(
    if 1:
        #table_html = ""
        #table_html += "<table border=\"1\">\n"
        for contig in rdata:
            page.tr()
            #table_html += "<tr>\n"
            index = 0
            for item in contig.split("\t"):
                index +=1
                if index == 8 or index == 7:
                    continue
                elif index == 1 and "File" not in item:
                    page.td("<a target=\"_blank\" href=\""+ref_asm+"\">%s</a>"%(ref_asm.split("/")[-1]))
                    continue

                page.td("<p style=\"font-size:15px\">"+item+"</p>")
                #table_html += "<td> %s </td> \n"%(item)
            page.tr.close()
            #table_html += "</tr>\n"
        page.table.close()
        #table_html += "</table>\n"
        #page.p( table_html )    
    #page.p.close()
    page.div.close()
    page.div( id_="quick", style_="background-color:%s;text-align:left;width:11%%;height:88%%;float:right;border:1px %s"%(backgroundColor, textColor))
    
    #items = ["<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>", ]
    #page.ul()
    page.p("<u>Quick summary</u><br>")
    page.p("#Reads:<br>%s"%readcount)
    page.p("#Contigs:<br>%s"%contigcount)
    page.p("#Scaffolds:<br>%s"%scaffoldcount)
    page.p("#ORFs:<br>%s"%orfcount)
    page.p("#Motifs:<br>%s"%motifcount)
    page.p("<a href=\"pipeline.summary\">Pipeline summary</a>")
    page.p("<a href=\"pipeline.commands\">Run commands</a>")
    page.p("<a href=\"https://github.com/treangen/metAMOS/wiki\">MetAMOS website</a>")
    page.div.close()
    page.div(id_="sideplots", style_="background-color:#FFFFFF;width:20.5%%;height:88%%;float:right;border:1px %s"%(textColor))
    #page.frameset(rows_="20%,20%,20%" ,cols_="100%")
    #page.div(id_="sideplot1", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    page.img(src_="ContigSizes.png",height_="25%",width_="100%")
    #page.div.close()
    #page.div(id_="sideplot2", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    page.img(src_="hist_contigs.png",width_="100%",height_="25%")
    #page.div.close()
    #page.div(id_="sideplot3", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    page.img(src_="ScaffoldSizes.png",width_="100%",height_="25%")
    #page.div.close()
    #page.div(id_="sideplot4", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    page.img(src_="hist_scaffold.png",width_="100%",height_="25%")
    page.div.close()
    #page.frameset.close()    
    #page.div.close()


    page.div( id_="krona", style="float:left;width:58%;height:74%")
    #page.iframe(src_="bar-phylum.png",style_="width:100%;height:100%;hspace=10")
    #page.frameset(rows_="100%" ,cols_="100%")
    page.div.close()

    #page.div( id_="metaphyler", style="float:left;width:28%;height:70%")
    #page.img(  hspace=10, alt='Abundance', src='bar-phylum.png' )
    #page.div.close()

    page.div(id="footer", style="background-color:#B8B8B8;clear:both;text-align:center;width:100%%;height:5%%;border:1px %s;vertical-align:middle"%(textColor))
    page.p("Generated %s"%(ds))
    page.div.close()
    #page.img( hspace=10, width=600, height=500, alt='Abundance', src='bar-phylum.png' )
    page.div.close()

    fout = open(html_prefix+"summary.html",'w')
    fout.write(page.__str__())
    fout.close()
