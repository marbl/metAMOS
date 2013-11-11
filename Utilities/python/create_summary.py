import sys, os, string, locale
ROOT = os.path.dirname(os.path.abspath(__file__))
#sys.path.insert(0, os.path.join(ROOT, '..'))
#sys.path.append(ROOT+"/lib")
import  markup
from pygooglechart import StackedVerticalBarChart
from pygooglechart import PieChart2D
from pygooglechart import PieChart3D
from pygooglechart import Axis
from pygooglechart import StackedHorizontalBarChart, StackedVerticalBarChart, \
         GroupedHorizontalBarChart, GroupedVerticalBarChart



import settings
from utils import *
import helper
from create_plots import *
from get_classify_stats import *
#let system set locale from available ones
_settings = Settings()

from datetime import datetime, date, time
locale.setlocale(locale.LC_ALL, '')

def intOrZero(string):
    if string:
        return locale.format("%d", int(string), grouping=True)
    else:
        return '0'

def getTable(header, content):
  result = "<div class=\"datagrid\"><table><thead><tr>"
  for (h, align) in header:
     result = "%s<th align=\"%s\">%s</th>"%(result, align, h)
  result = "%s</tr></thead><tbody>"%(result)
  
  i = 0
  for row in content:
     if i % 2 == 0:
        result = "%s<tr>"%(result)
     else:
        result = "%s<tr class=\"alt\">"%(result)
     for (col, align) in row:
        result = "%s<td align=\"%s\">%s</td>"%(result, align, col)
     result = "%s</tr>"%(result)
     i += 1
  result = "%s</tbody></table></div>"%(result)
  return result

def outputLibraryInfo(headerArray, dataArray, outputHeader, libcnt, format, mated, interleaved, mmin, mmax, outputFastQC):
   if outputHeader:
      headerArray.append(["Library #", "left"])
      headerArray.append(["Format", "left"])
      headerArray.append(["Mated?", "left"])
      headerArray.append(["Min Insert Size", "right"])
      headerArray.append(["Max Insert Size", "right"])
      if outputFastQC:
         headerArray.append(["First FastQC Report", "left"])
         headerArray.append(["Second FastQC Report", "left"]) 

   row = []
   row.append([libcnt, "left"])
   row.append([format, "left"])
   row.append([mated, "left"])
   row.append([mmin, "right"])
   row.append([mmax, "right"])

   if outputFastQC and os.path.exists("lib%d.1.fastqc/fastqc_report.html"%(libcnt)):
      if mated.lower() == "true":
         if interleaved.lower() == "true":
            row.append(["<a target=\"_blank\" href=\"lib%d.1.fastqc/fastqc_report.html\">interleaved</a>"%(libcnt), "left"])
            row.append(["NA", "left"])
         else:
            row.append(["<a target=\"_blank\" href=\"lib%d.1.fastqc/fastqc_report.html\">left</a>"%(libcnt), "left"])
            row.append(["<a target=\"_blank\" href=\"lib%d.2.fastqc/fastqc_report.html\">right</a>"%(libcnt), "left"])
      else:
         row.append(["<a target=\"_blank\" href=\"lib%d.1.fastqc/fastqc_report.html\">unmated</a>"%(libcnt), "left"])
         row.append(["NA", "left"])
   dataArray.append(row)

def outputValidate(headerArray, dataArray, outputHeader, best, results):
   assembler = results[0]
   isFirst = True
   if outputHeader:
      for r in results:
         if isFirst:
            headerArray.append([r.title(), "left"])
            isFirst = False
         else:
            headerArray.append([r.title(), "right"])
      return

   isBest = False
   if best.lower() == assembler.lower():
      isBest = True

   row = []
   counter = 0
   isFirst = True
   for r in results:
      if isFirst:
         row.append(["%s%s%s"%("<b>" if isBest else "", r, "</b>" if isBest else ""), "left"])
         isFirst = False
      else:
         score = ""
         if r == None or r.lower() == "none":
            score = "N/A"
         else:
            if headerArray[counter][0].upper() in SCORE_TYPE.mapping:
               scoreNum = SCORE_TYPE.mapping[headerArray[counter][0].upper()]
               if scoreNum == SCORE_TYPE.LAP or scoreNum == SCORE_TYPE.ALE or scoreNum == SCORE_TYPE.CGAL:
                  score = "%.6f"%(float(r)) 
               elif scoreNum == SCORE_TYPE.SNP or scoreNum == SCORE_TYPE.FRCBAM or scoreNum == SCORE_TYPE.ORF:
                  score = "%d"%(int(r))
               elif scoreNum == SCORE_TYPE.REAPR:
                  score = "%.2f"%(float(r))
               elif scoreNum == SCORE_TYPE.N50:
                  score = "{:,}".format(int(r))
            elif len(r) != 0:
               score = "%.4f"%(float(r))
         row.append(["%s%s%s"%("<b>" if isBest else "", score, "</b>" if isBest else ""), "right"])
      counter += 1
   dataArray.append(row)

def create_summary(first,amosbnk,prefix,ref_asm,utils,img,rund,nLibs,taxa_level,dbdir):
#if __name__ == "__main__":
#    if len(sys.argv) < 5:
#        print "usage: create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm> <Utils dir> <run dir> <# of libs> <taxa level of classifications>"
#        sys.exit(0)
    
    html_prefix = prefix
    prefix = prefix.replace("/html/", "")
    MA_dir = prefix
    MA_dir = MA_dir.replace("/Postprocess/out","")
    mp = open(first,'r')
    #mp2 = open(sys.argv[1].replace("s12","s3"),'r')    

    # set working dir
    os.chdir(html_prefix)

    if not os.path.exists(html_prefix+"asmstats.out"):
        libPath = rund.replace("bin", "lib")
        #print "perl -I %s %s/perl/statistics.pl %s > %sasmstats.out"%(libPath,utils,ref_asm,html_prefix)
        run_process(_settings,"perl -I %s %s/perl/statistics.pl %s > %sasmstats.out"%(libPath,utils,ref_asm,html_prefix),"Classify")
    report = open(html_prefix+"asmstats.out",'r')

    initialStep = "Annotate"

    # get metamos version
    version = "1.0"
    summary = open("%s/pipeline.run"%(MA_dir), 'r')
    for line in summary:
       line = line.replace("\n","")
       if "#" in line:
          continue
       elif "metAMOS Version:" in line:
          version = line.replace("\n","").split("\t")[-1]
    summary.close()

    steps = []
    steps.append("Preprocess")
    steps.append("Assemble")
    steps.append("MapReads")
    steps.append("Validate")
    #steps.append("MultiAlign")
    steps.append("FindORFS")
    steps.append("FindRepeats")
    steps.append("Scaffold")
    steps.append("FindScaffoldORFS")
    steps.append("Abundance")
    steps.append("Annotate")
    steps.append("FunctionalAnnotation")
    steps.append("Propagate")
    steps.append("Classify")
    steps.append("Browse Results")

    step_status = {}
    step_status["Preprocess"] = "OK"
    step_status["Assemble"] = "OK"
    step_status["MapReads"] = "OK"
    step_status["Validate"] = "OK"
    #step_status["MultiAlign"] = "OK"
    step_status["FindORFS"] = "OK"
    step_status["FindRepeats"] = "OK"
    step_status["Scaffold"] = "OK"
    step_status["FindScaffoldORFS"] = "OK"
    step_status["Abundance"] = "FAIL"
    step_status["Annotate"] = "FAIL" 
    step_status["FunctionalAnnotation"] = "FAIL" 
    step_status["Propagate"] = "NONE"
    step_status["Classify"] = "NONE"
    step_status["Browse Results"] = "OK"
    ##get status of each step from Log dir
    for step in steps:
        if step == "Browse Results":
            continue
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
        if completed and not skipped:
            step_status[step] = "OK"
        elif skipped:
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
    #os.system("python %s/python/create_plots.py %s/plot.tab proba1"%(utils,html_prefix))
    create_plots("%s/plot.tab"%(html_prefix),"%s"%("proba1"))

    ##update counts
    #count reads
    os.system("cat `ls %s/Preprocess/out/* |egrep 'lib[0-9]*\.fasta$'` | grep -c \">\" > readcount.txt"%(MA_dir))
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
    os.system("grep -c \">\" %s/FindORFS/out/proba.fna > orfcount.txt"%(MA_dir))
    orfcount = open("orfcount.txt",'r').read().replace("\n","")  
    #print orfcount
    os.system("rm orfcount.txt")
    ##copy stuff
    for step in steps:
#        step = step.lower()
        if os.path.exists("%s/javascript/%s.html"%(utils,step)):
            os.system("cp %s/javascript/%s.html %s/"%(utils,step,html_prefix))
    os.system("cp %s/Logs/COMMANDS.log %s/pipeline.commands"%(MA_dir,html_prefix))
    os.system("cp %s/pipeline.run %s/pipeline.summary"%(MA_dir,html_prefix))
    os.system("cp %s/javascript/style.css %s/"%(utils,html_prefix)) # TEMP: change back to cp
    os.system("cp -r %s/../KronaTools/src %s/../KronaTools/img %s/"%(utils, utils, html_prefix)) # TODO: unhack KronaTools path
#    os.system("cp %s/blocks.jpg %s/"%(img,html_prefix))
#    os.system("cp %s/blocks_small.jpg %s/"%(img,html_prefix))
    os.system("cp %s/blocks_dark_tiny.png %s/"%(img,html_prefix))
    os.system("cp %s/name.png %s/"%(img,html_prefix))

    # generate dynamic java scripts
    # first classify and propagate
    #os.system("python %s/python/get_classify_stats.py %s/propagate.in.clusters %s/propagate.out.clusters %s/DB/tax_key.tab %s Classify.html Propagate.html %s"%(utils, html_prefix, html_prefix, utils, html_prefix, taxa_level)) 
    get_classify_stats("%s/propagate.in.clusters"%(html_prefix),"%s/propagate.out.clusters"%(html_prefix),"%s/tax_key.tab"%(dbdir),"%s"%(html_prefix),"Classify.html","Propagate.html","%s"%(taxa_level))

    # generate preprocess
    preprocess = markup.page()
    preprocess.init(css="style.css")#bodyattrs={'style':"margin:0px"})
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
    headerArray = []
    dataArray = []
    for line in summary:
       line = line.replace("\n","")
       if "#" in line:
          continue
       elif "asmcontigs:" in line:
          asmc = line.replace("\n","").split("\t")[-1]
          if len(asmc) <= 2:
             continue
          preprocess.add("<div class=\"datagrid\">Pre-assembled contigs input: %s</div>"%(asmc))
          preprocess.br()
       elif "format:" in line:
          if format and not libadded:
             outputLibraryInfo(headerArray, dataArray, firstLib, libcnt, format, mated, interleaved, mmin, mmax, nQC > 0)
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

       elif "frg" in line:
          data = line.split("\t")
    if format and not libadded:
       outputLibraryInfo(headerArray, dataArray, firstLib, libcnt, format, mated, interleaved, mmin, mmax, nQC > 0)
    preprocess.add(getTable(headerArray, dataArray))
    summary.close()

    if os.path.exists("%s/Postprocess/out/kmergenie_report.html"%(MA_dir)):
       preprocess.iframe(id_="KmerGenie", src_="%s/Postprocess/out/kmergenie_report.html"%(MA_dir), width="800", height="5000")

    preprocess_out = open("%s/Preprocess.html"%(html_prefix), 'w')
    preprocess_out.write(preprocess.__str__())
    preprocess_out.close()

    validate_out = open("%s/Validate.html"%(html_prefix), 'w')
    headerArray = []
    dataArray = []
    firstScore = True
    bestAsm = ""
    refs = [] 
    if os.path.exists("%s/Postprocess/out/lap.scores"%(MA_dir)):
       best = open("%s/Postprocess/out/best.asm"%(MA_dir), 'r')
       bestAsm = best.read()
       best.close()
       laps = open("%s/Postprocess/out/lap.scores"%(MA_dir), 'r')
       validate = markup.page()
       validate.init(css="style.css")
       validate.p()
       validate.add("<div class=\"datagrid\">")
       validate.add("Selected assembler: %s"%(bestAsm))
       ref = open("%s/Postprocess/out/ref.asm"%(MA_dir), 'r')
       first = True
       for r in ref.xreadlines():
           refs.append(r)
           if first:
              validate.br()
              first = False
           validate.add("Selected reference: %s"%(r))
       ref.close()
       if os.path.exists("%s/Classify/out/contaminant.true"%(MA_dir)):
          cont = open("%s/Classify/out/contaminant.true"%(MA_dir), 'r')
          contPercent = cont.read().split()
          cont.close()
          validate.br()
          validate.add("<b>Sample may have contaminants, only %s%% assigned to %s. Check Annotate output</b>"%(contPercent[0], contPercent[1]))
       validate.br()
       for line in laps:
          line = line.replace("\n","")
          if "#" in line:
             continue
          else:
             res = line.split("\t")
             outputValidate(headerArray, dataArray, firstScore, bestAsm, res)
             firstScore = False
       validate.add(getTable(headerArray, dataArray))
       validate.add("</div>")
       laps.close()

    # when we have quast, we will add our table to their report, otherwise write standalone report
    if os.path.exists("%s/Postprocess/out/quast/report.html"%(MA_dir)):
       os.system("cp %s/Postprocess/out/quast/report.html %s/Postprocess/out/quast/~report.html"%(MA_dir, MA_dir))
       quastIn = open("%s/Postprocess/out/quast/~report.html"%(MA_dir), 'r')
       quastOut = open("%s/Postprocess/out/quast/report.html"%(MA_dir), 'w')
       skip = False
       for line in quastIn.xreadlines():
          if not skip:
             quastOut.write(line + "\n")
          else:
             quastOut.write("margin-left: 10px;")
             skip = False

          if ".content" in line:
             skip = True
       quastIn.close()
       quastOut.close()
       validate.iframe(id_="quast", src_="%s/Postprocess/out/quast/report.html"%(MA_dir), width="800", height="1000")
    validate_out.write(validate.__str__())
    validate_out.close()

    # multialign step
    treeScripts= {} 
    treeScripts["http://www.jsphylosvg.com/js/jquery/jquery-1.4.2.min.js"] = "javascript"
    treeScripts["http://www.jsphylosvg.com/js/jquery/jquery.simplemodal.1.4.1.min.js"] = "javascript"
    treeScripts["http://www.jsphylosvg.com/js/raphael/raphael-min.js"] = "javascript"
    treeScripts["http://www.jsphylosvg.com/js/jsphylosvg-min.js?1.29"] = "javascript"
    treeScripts["http://www.jsphylosvg.com/js/unitip/js/unitip.js"] = "javascript"

    treeScript = [] 
    treeScript.append("<script type=\"text/javascript\">")
    treeScript.append("   var dataObject = { newick: '(M_canettii_CIPT_d:0.059437,(M_canettii_CIPT_e:0.07845,((((M_africanum_GM04118:0.00646,(M_bovis_AF2122_97:0.00269,(M_bovis_BCG_Tokyo_1:0.00013,((M_bovis_BCG_Korea_1:0.0,M_bovis_BCG_Pasteur:0.0):0.00014,M_bovis_BCG_Mexico:0.00013):0.00014):0.00294):0.00514):0.00207,(((MTBC_Beijing_NITR20:0.00914,(MTBC_CCDC5180:0.00102,(MTBC_CCDC5079:0.00423,MTBC_CCDC5079:0.00014):0.00078):0.00016):0.00368,((((MTBC_Erdman_ATCC_35:0.00075,MTBC:0.00123):0.00129,(MTBC_CDC1551:0.00245,MTBC_Haarlem3_NITR2:0.01398):0.00093):0.00065,((((MTBC_CTRI_2:0.001,(MTBC_KZN_4207:0.00014,(MTBC_KZN_605:0.00012,MTBC_KZN_1435:0.00014):0.00012):0.00071):0.00061,(MTBC_RGTB327:0.00462,MTBC_UT205:0.0012):0.00037):0.00014,MTBC_F11:0.00158):0.00109,(MTBC_H37Ra:0.00016,(MTBC_H37Rv:0.00028,MTBC_H37Rv:0.00014):0.00014):0.00281):0.00047):0.00112,MTBC_CAS_NITR204:0.009827):0.00022):0.00247,((MTBC_EAI5:0.0023,MTBC_EAI5_NITR206:0.00933):0.0035,(MTBC_RGTB423:0.00868,spades.45.asm.contig:0.0035):0.00066):0.00289):0.00015):0.07224,(M_canettii_CIPT_a:0.01891,M_canettii_CIPT_b:0.05096):0.03720):0.03909,M_canettii_CIPT_c:0.056869):0.00469):0.02099);' };")
    treeScript.append("   function load() {")
    treeScript.append("      var divVal = $('#svgCanvas').empty();")
    treeScript.append("      var isCircular = $('#circularize-value')[0].value;")
    treeScript.append("      Smits.PhyloCanvas.Render.Parameters.Rectangular.bufferX = 150;")
    treeScript.append("      Smits.PhyloCanvas.Render.Parameters.Circular.bufferRadius = 0.35;")
    treeScript.append("      Smits.PhyloCanvas.Render.Style.line.stroke = 'rgb(0,0,255)';")
    treeScript.append("      Smits.PhyloCanvas.Render.Style.text[\"font-size\"] = 10;")
    treeScript.append("      if (isCircular.toLowerCase() == 'true') {")
    treeScript.append("         phylocanvas = new Smits.PhyloCanvas(")
    treeScript.append("                                         dataObject,")
    treeScript.append("                                         'svgCanvas',")
    treeScript.append("                                          1000, 1000, 'circular'")
    treeScript.append("                                         );")
    treeScript.append("         $('#circularize-value')[0].value = 'false';")
    treeScript.append("      } else {")
    treeScript.append("         phylocanvas = new Smits.PhyloCanvas(")
    treeScript.append("                                         dataObject,")
    treeScript.append("                                         'svgCanvas',")
    treeScript.append("                                          500, 500")
    treeScript.append("                                         );")
    treeScript.append("         $('#circularize-value')[0].value = 'true';")
    treeScript.append("      }")
    treeScript.append("      init();")
    treeScript.append("   }")
    treeScript.append("function snapshot() {")
    treeScript.append("      var svgSource = phylocanvas.getSvgSource();")
    treeScript.append("      if(svgSource) {")
    treeScript.append("         var url = 'data:image/svg+xml,' + svgSource")
    treeScript.append("         window.open(url);")
    treeScript.append("      }")
    treeScript.append("   }")
    treeScript.append("</script>")
    treeBody = {}
    treeBody["onload"] = "load()"

    maStep = markup.page()
    maStep.add("\n".join(treeScript))
    maStep.init(
                css=("http://www.jsphylosvg.com//js/yui/build/cssfonts/fonts-min.css", "http://www.jsphylosvg.com/js/unitip/css/unitip.css"),
                script=treeScripts, 
                bodyattrs=treeBody)
    maStep.div()
    maStep.add("<input type=\"button\" id=\"download-link\" onclick=\"snapshot()\" value=\"Snapshot\" />")
    maStep.div.close()
    maStep.div()
    maStep.add("<input type=\"button\" id=\"circularize-link\" onclick=\"load()\" value=\"Toggle Circular\" />")
    maStep.add("<input type=\"hidden\" id=\"circularize-value\" value=\"false\">")
    maStep.div.close()
    maStep.div(id_="svgCanvas")
    maStep.div.close()
    maStep_out = open("%s/MultiAlign.html"%(html_prefix), 'w')
    maStep_out.write(maStep.__str__())
    maStep_out.close()

    # todo, need to add report for MapReads including # reads mapped (%), contig coverage histogram, and % reads between contigs and number of links histogram. Also re-estimated insert sizes for each lib
    #mapreads = markup.page()
    #mapreads.init(bodyattrs={'style':"margin:0px"})
    #mapreads.img(src_="hist_ctgcvg.png",height_="100%",width_="100%")
    #mapreads_out = open("%s/MapReads.html"%(html_prefix), 'w')
    #mapreads_out.write(mapreads.__str__())

    ##This will create ScaffoldSizes.png,ContigSizes.png
    
    ##create code to automatically generate .js files for HTML report
    ## let's start with Abundance
    #abundance_js = open(html_prefix+"abundance.js",'w')
    
    rdata = []
    for line in report:
        rdata.append(line)
       
    if not os.path.exists(html_prefix+"covstats.out"):
        run_process(_settings,"%s/analyze-read-depth -x 2 %s > %scovstats.out"%(rund,amosbnk,html_prefix),"Classify")
    ff = open(html_prefix+"covstats.out",'r')
    covdata = []
    #covdata = ff.readlines()
    #zflag = 0
    for line in ff:
        covdata.append(line)
    
    if not os.path.exists(html_prefix+"stats.out"):
        #os.system("%s/astats %s > %sstats.out"%(rund,amosbnk,html_prefix))
        run_process(_settings,"%s/astats %s > %sstats.out"%(rund,amosbnk,html_prefix),"Classify")
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
    ds = dt.strftime("%d %b %Y,<br/>%l:%M%p")
    title = "metAMOS report"
    # header is blank for now
    header = []

    # write the javascript we need on the page
    script = []
    script.append("<script type=\"text/javascript\">")
    script.append("var steps = ['%s'];"%("','".join(steps)))
    script.append("function load(step) {")
    for step in steps:
       if step_status[step] != "OK":
          script.append("   if (step.toLowerCase() == \"%s\".toLowerCase()) {"%(step))
          script.append("      return;")
          script.append("   }")
    script.append("   for (var i = 0; i < steps.length; i++) {")
    script.append("      var current = steps[i].toLowerCase();")
    script.append("      if (current == step) {")
    script.append("         document.getElementById(step + 'Button').className = 'menuItemSelected';")
    script.append("         document.getElementById(step).style.display = 'block';")
    script.append("         document.getElementById(step + 'Marker').className = 'markerSelected';")
    script.append("      } else {")
    script.append("         document.getElementById(current + 'Button').className = 'menuItemUnselected';")
    script.append("         document.getElementById(current).style.display = 'none';")
    script.append("         document.getElementById(current + 'Marker').className = 'markerUnselected';")
    script.append("      }")
    script.append("   }")
    script.append("}")
    script.append("")
    script.append("</script>")

    #<link rel="shortcut icon" href="../assets/ico/favicon.ico">
    #<link rel="apple-touch-icon-precomposed" sizes="144x144" href="../assets/ico/apple-touch-icon-144-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" sizes="114x114" href="../assets/ico/apple-touch-icon-114-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" sizes="72x72" href="../assets/ico/apple-touch-icon-72-precomposed.png">
    #<link rel="apple-touch-icon-precomposed" href="../assets/ico/apple-touch-icon-57-precomposed.png">

    # generate the dictionary of javascript pages we need
    scripts = {}
#    scripts["file://%s/javascript/jquery-latest.js"%(utils)] = "javascript"
#    scripts["http://code.jquery.com/jquery-latest.js"] = "javascript"
    # now a javascript for each page
    for step in steps:
       scripts["%s.js"%(step.lower())] = "javascript"

    # attributes for the body tax
    body = {}
    body["onload"] = "load('%s')"%(initialStep.lower())

    # the footer for the page

    footer = ""#"Generated %s"%(ds)
    #styles = ( 'style2.css')#'./html/bootstrap.css', './html/boostrap-responsive.css')#'style2.css')#'layout.css', 'alt.css', 'images.css' )
    #styles = ( 'layout.css', 'alt.css', 'images.css' )
    #meta = ('viewport':"width=device-width, initial-scale=1.0",'description':'','author':'')

    # create HTML page
    page = markup.page( )
    # add the javascript free form
    page.add("\n".join(script))
    # initialize the rest of the body/html headers
    page.init( title=title,   \
               css='style.css', \
               script=scripts,\
               header='\n'.join(header), \
               bodyattrs=body, \
               footer=footer )
    
#    page.br()
    #page.div( class_ = 'navbar navbar-fixed-top')
    #page.div( class_ = 'navbar-inner')
    
    #page.div( class_ = 'container', style_='width:100%;height=100%')
    
    #page.div.close()
    #page.div.close()
    #page.div.close()
    #page.div( id_='page' )
    
    # header
    
#    page.add('<table style="width:100%;"><tr><td><img src="blocks_small.jpg"/></td><td><a target="_blank" href="https://github.com/treangen/metAMOS/wiki"><div><img src="name.jpg"/><br/>Under peer review</a></div></td><td class="title" style="width:100%"></td></tr></table>')
    
    '''
    page.div( id_='header', style="background-color:#B8B8B8;clear:both;text-align:center;width:100%;height:7%;border:1px solid black") 
    page.h1("<u>MetAMOS <font color=\"blue\">v1.0</font> metagenomic assembly & analysis report</u>" , style_="text-align:center;vertical-align:top")
    page.h1.close()
    #page.div( id_='title')# style="vertical-align:bottom;")
    #page.font( size=14)
    #page.br( )
    page.div.close()
    '''
    
    #<frameset rows="40%,60%"cols="80%,20%">
    #<frame src="report.krona.html">
    #<frameset rows="20%,20%">
    #<frame src="report.krona.html">
    #<frame src="report.krona.html">
    #</frameset>
    #<frameset rows="60%" cols="40%,40%">
    #<frame src="report.krona.html">
    #<frame src="report.krona2.html">
    #</frameset>f
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
    
    page.table(style_="width:100%;height:100%;")
    page.tr()
    page.td( id_="menu", style_="padding:0px;")
    page.div(style_="height:100%;")
    
#    page.table(style_="height:100%")
#    page.tr()
#    page.td()#style_="border-right:1px solid black")
#    page.div(style_="box-shadow:inset -1px -1px 5px #555555;")
    page.add("<a target=\"_blank\" href=\"https://github.com/treangen/metAMOS/wiki\"><img style=\"padding-top:5px;\" src=\"name.png\"/>")
    page.add("<img src=\"blocks_dark_tiny.png\"/></a>")
#    page.add("<img src=\"blocks_tiny2.jpg\"/>")
    page.add("<div style=\"padding:2px;font-size:12px;\"><a target=\"_blank\" href=\"http://genomebiology.com/2013/14/1/R2\">Treangen TJ, Koren S, et. al.  Genome Biol. 2013 Jan 15;14(1):R2. PMID: 23320958.</a></div>")
    page.add("<br/>")
#    page.div.close()
#    page.td.close()
#    page.tr.close()
    
    #items = ["<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>", ]
    #page.ul()
#    page.tr()
#    page.td(style_="padding:0px;")
    page.table( id_="links", class_="menu" )
    for step in steps:
       page.tr()
       status = "NA"
       try:
          if step_status[step] == "OK":
             status = "OK"
          elif step_status[step] == "FAIL":
             status = "FAIL"
          elif step_status[step] == "SKIP":
             status = "SKIPPED"
       except KeyError:
          continue

       tableHTML = []
       tableHTML.append("<td class=\"menuItemTd\"><a class=\"menuLink\" href=\"javascript:void(0)\"><div class=\"menuItem %s\"><div id=\"%sButton\" class=\"menuItemUnselected\" onclick=\"load('%s');\">"%(status.lower(), step.lower(), step.lower()))
       tableHTML.append("<table style='width:100%%;'><tr><td style='width:100%%;'><span class='step %sText'>%s</span>"%(status.lower(), step))
       tableHTML.append("<br>")
       tableHTML.append("<span class='status %sText'>%s</span></td></tr></table>"%(status.lower(), status))
       tableHTML.append("</div></div></a></td><td><div id='%sMarker' class='markerUnselected'>&#x25CF;</div></td>"%(step.lower()))
       page.add("\n".join(tableHTML))
       page.tr.close()

    page.table.close()
#    page.td.close()
#    page.tr.close()
    #page.ul.close()
    
#    page.tr()
#    page.td(style_="height:100%;")#border-right:1px solid black")
    page.div(class_="notes")
    page.add("<div style=\"font-size:14px;\"><br/><a target=\"_blank\" href=\"pipeline.summary\">Pipeline summary</a><br/>")
    page.add("<a target=\"_blank\" href=\"pipeline.commands\">Run commands</a><br/><br/></div>")
    tableHTML = []
    tableHTML.append("<table style=\"font-size:12px\"><tr>")
    tableHTML.append("<tr><td>Version:</td><td>%s</td></tr>"%(version))
    tableHTML.append("<tr><td>Created:<br/>&nbsp;</td><td>%s</td></tr>"%(ds))
    tableHTML.append("</table>")
    page.add("\n".join(tableHTML))
    page.div.close()
#    page.td.close()
#    page.tr.close()
#    page.table.close()
    
    page.div.close()
    page.td.close()
    
    page.td( id_="krona", style_="width:100%;height:100%padding:0px;")
    page.table(class_="charts",style_="width:100%;height:100%")
    page.tr()
    page.td(class_="main", style_="height:100%", colspan_="4")
    page.div(class_="shadow")
    page.table(style_="width:100%;height:100%")
    page.tr()
    page.td(class_="corner")
    page.td.close()
    page.td()
    page.td.close()
    page.td(class_="corner")
    page.td.close()
    page.tr.close()
    page.tr()
    page.td()
    page.td.close()
    page.td(class_="inset")
    for step in steps:
        if step == "Browse Results":
            page.iframe( id_=step.lower(), src_="../.", style_="display:none;" )
        else:
            page.iframe( id_=step.lower(), src_="%s.html"%(step), style_="display:none;" )
        page.iframe.close()
    page.td.close()
    page.td()
    page.td.close()
    page.tr.close()
    page.tr()
    page.td(class_="corner")
    page.td.close()
    page.td()
    page.td.close()
    page.td(class_="corner")
    page.td.close()
    page.tr.close()
    page.table.close()
    page.div.close()
    page.td.close()
    page.tr.close()
    #page.tr()
    #page.div(id_="sideplots")#, style_="background-color:#FFFFFF;width:20.5%%;height:88%%;float:right;border:1px")
    #page.div(id_="sideplot1", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    #page.td(style_="width:25%")
    #page.img(class_="chart", src_="ContigSizes.png")
    #page.td.close()
    #page.td(style_="width:25%")
    #page.div.close()
    #page.div(id_="sideplot2", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    #page.img(class_="chart", src_="hist_contigs.png")
    #page.td.close()
    #page.td(style_="width:25%")
    #page.div.close()
    #page.div(id_="sideplot3", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    #page.img(class_="chart", src_="ScaffoldSizes.png")
    #page.td.close()
    #page.td(style_="width:25%")
    #page.div.close()
    #page.div(id_="sideplot4", style_="background-color:#FFFFFF;width:20.5%;height:22%;float:right")
    #page.img(class_="chart", src_="hist_scaffold.png")
    #page.td.close()
    #page.div.close()
    #page.div.close()
    #page.td.close()
    #page.tr.close()
    page.table.close()
    #page.iframe(src_="bar-phylum.png",style_="width:100%;height:100%;hspace=10")
    #page.frameset(rows_="100%" ,cols_="100%")
    page.td.close()
    
    #page.li("<a target=\"_blank\" href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>")    
    #page.li("<a href=\"http://cbcb.umd.edu/~mpop/Software.shtml\">Related software</a>")
    #page.li("<a href=\"http://cbcb.umd.edu\">CBCB</a>")
    
#    page.li( items[1] )
    #page.ul.close( )
    #page.div.close()
    
    page.div( id_="content")

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
    '''
    #page.table.close()        
    #page.p( paragraphs )
    #page.p(style_="font-size:6px")
    page.table(border="0")
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
    '''
    page.td( id_="quick" )
    page.div(style_="height:100%")
    #items = ["<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>", ]
    #page.ul()
    page.br()
    page.div(class_="stats")
    page.table(class_="stats")
    page.add("<tr><td class=\"number\">%s</td><td>Reads</td></tr>"%intOrZero(readcount))
    page.add("<tr><td class=\"number\">%s</td><td>Contigs</td></tr>"%intOrZero(contigcount))
    page.add("<tr><td class=\"number\">%s</td><td>Scaffolds</td></tr>"%intOrZero(scaffoldcount))
    page.add("<tr><td class=\"number\">%s</td><td>ORFs</td></tr>"%intOrZero(orfcount))
    page.add("<tr><td class=\"number\">%s</td><td>Motifs</td></tr>"%intOrZero(motifcount))
    page.table.close()
    page.div.close()

    page.div.close()
    page.td.close()
    page.tr.close()
    page.table.close()

    #page.div( id_="metaphyler", style="float:left;width:28%;height:70%")
    #page.img(  hspace=10, alt='Abundance', src='bar-phylum.png' )
    #page.div.close()

    #page.img( hspace=10, width=600, height=500, alt='Abundance', src='bar-phylum.png' )
    
    #page.div.close()

    fout = open(html_prefix+"summary.html",'w')
    fout.write(page.__str__())
    fout.close()
