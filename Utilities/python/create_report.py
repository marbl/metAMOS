import sys, os, string
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, '..'))
sys.path.append(ROOT+"/lib")
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

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "usage: create_report.py <metaphyler tab file> <AMOS bnk> <output prefix> <ref_asm>"
        sys.exit(0)
    ref_asm = sys.argv[4]
    mp = open(sys.argv[1],'r')
    mp2 = open(sys.argv[1].replace("s12","s3"),'r')    
    if not os.path.exists("./out/asmstats.out"):
        os.system("perl ./lib/statistics.pl %s > ./out/asmstats.out"%(sys.argv[4]))        
    report = open("./out/asmstats.out",'r')
    
    rdata = []
    for line in report:
        rdata.append(line)
    #if not os.path.exists("./out/asmstats2.out"):
    #    os.system("perl ./lib/statistics.pl ./even/s3.fa > ./out/asmstats2.out")#%(sys.argv[5]))        
    #report = open("./out/asmstats2.out",'r')
    #report.readline()
    #for line in report:
    #    rdata.append(line)        
       
    if not os.path.exists("./out/covstats.out"):
        os.system("analyze-read-depth -i %s -x 10 -l 500 > ./out/covstats.out"%(sys.argv[2]))
    ff = open("./out/covstats.out",'r')
    covdata = []
    #covdata = ff.readlines()
    #zflag = 0
    for line in ff:
        covdata.append(line)
    
    if not os.path.exists("./out/stats.out"):
        os.system("astats %s > ./out/stats.out"%(sys.argv[2]))
    dd = open("./out/stats.out",'r')
    ddata = dd.readlines()
    prefix = sys.argv[3]
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
    for line in mpd:
        data = line.split("\t")
         
        ids.append(data[0])
        cov.append(int(float(data[2])))
        abund.append(float(data[1]))

    mp2.readline()
    cov2 = []
    abund2 = []
    ids2 = []
    mpd = mp2.readlines()
    for line in mpd:
        data = line.split("\t")
        ids2.append(data[0])
        cov2.append(int(float(data[2])))
        abund2.append(float(data[1]))        
        
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
    chart2.download('./out/abund.png')

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
    chart.download('./out/bar-phylum.png')
    
    dt = datetime.now()
    ds = dt.strftime("%A, %d. %B %Y %I:%M%p")
    title = "metAMOS: a metagenomic assembly pipeline for AMOS"
    header = ""#metAMOS Metagenomic assembly report"

    footer = ""#"Generated %s"%(ds)
    styles = ( 'style2.css')#'layout.css', 'alt.css', 'images.css' )
    #styles = ( 'layout.css', 'alt.css', 'images.css' )

    page = markup.page( )
    page.init( css=styles, title=title, header=header, footer=footer )

#    page.br()
    page.div( id_='page' )
    page.div( id_='header' )
    page.div( id_='title')# style="vertical-align:bottom;")
    #page.font( size=14)
    page.p("metAMOS metagenomic assembly report")

    #page.br( )
    page.div.close()
    
    #font, estimate
    page.div.close()
    page.div( id_="menu")
    items = ["<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>", ]
    page.ul()
    page.li("<a href=\"http://cbcb.umd.edu/software/metamos\">metAMOS website</a>")    
    page.li("<a href=\"http://cbcb.umd.edu/~mpop/Software.shtml\">Related software</a>")
    page.li("<a href=\"http://cbcb.umd.edu\">CBCB</a>")
#    page.li( items[1] )
    page.ul.close( )
    page.div.close()
    page.div( id_='wrapper')
    #page.div( id_="content")
    
    #page.a( "Reference assembly:", class_='internal', href='%s'%(ref_asm) )          

    #paragraphs = ( ddf )
    cnt = 0
    #table_html = ""
    if 0:
        page.table(border="1")
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

    page.table(border="1")
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
                    page.td("<a href=\""+ref_asm+"\">%s</a>"%(ref_asm.split("/")[-1]))
                    continue
                page.td(item)
                #table_html += "<td> %s </td> \n"%(item)
            page.tr.close()
            #table_html += "</tr>\n"
        page.table.close()
        #table_html += "</table>\n"
        #page.p( table_html )    

    page.img( hspace=10, width=600, height=500, alt='Abundance', src='bar-phylum.png' )
    page.p("Generated %s"%(ds))
#    page.p("Metaphyler predicted abundance by class")
#    page.img( hspace=50, width=600, height=300, alt='Abundance', src='class.png' )
 #   page.br()    
  #  page.p("Metaphyler predicted abundance by genus")
  #  page.img( hspace=50, width=600, height=300, alt='Abundance', src='cov.png' )
  #  page.br()    
    page.div.close()
    page.div.close()

    fout = open(prefix+".html",'w')
    fout.write(page.__str__())
    fout.close()
