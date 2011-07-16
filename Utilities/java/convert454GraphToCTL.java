import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;

/* convert this format
* 1       contig00001     694     5.7
* 2       contig00002     390     3.7
* ...
* C       25      3'      26      3'      4
* C       26      3'      33      3'      5
* C       32      3'      34      5'      3
* C       56      5'      2621    3'      2
* C       56      5'      4649    5'      3
* to AMOS CTE messages
{CTE
adj:N
nds:0,6
iid:2
sze:2512
lnk:
21377
6163
.
std:212
}
*/
class convert454GraphToCTL {
   public static void main(String[] args) throws Exception {
      boolean outputCTE = false;
      if (args.length > 2) { 
         outputCTE = Boolean.parseBoolean(args[2]);
      }

      HashMap<String, String> eidTOiid = new HashMap<String, String>();
      HashMap<String, String> idToUID = new HashMap<String, String>();
      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(args[0])));
      String line = null;
      while ((line = bf.readLine()) != null) {
         String[] splitLine = line.trim().split("\\s+");
          System.err.println("STORING ID " + splitLine[0] + " TO " + splitLine[1]);
         eidTOiid.put(splitLine[0], splitLine[1]);
      }
      bf.close();

      HashMap<String, String> ctes = new HashMap<String, String>();
      bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(args[1])));
      line = null;
      int iid = 1;
      while ((line = bf.readLine()) != null) {
         String[] splitLine = line.trim().split("\\s+");
         if (line.startsWith("C")) {
           String ctesID = splitLine[1] + "_" + splitLine[3];
           String orient = "";
           StringBuffer cte = new StringBuffer();
           cte.append("{CTE\n");
           if ("3'".equalsIgnoreCase(splitLine[2]) && "3'".equalsIgnoreCase(splitLine[4])) {
              orient = "I";
           } else if ("5'".equalsIgnoreCase(splitLine[2]) && "3'".equalsIgnoreCase(splitLine[4])) {
              orient = "A";
           } else if ("3'".equalsIgnoreCase(splitLine[2]) && "5'".equalsIgnoreCase(splitLine[4])) {
              orient = "N";
           } else if ("5'".equalsIgnoreCase(splitLine[2]) && "5'".equalsIgnoreCase(splitLine[4])) {
              orient = "O";
           } else {
              System.err.println("UNRECOGNIZED ORIENATATION " + line);
              System.exit(1);
           }
           ctesID = ctesID + "_" + orient;
           if (ctes.get(ctesID) != null) {
             System.err.println("ERROR ID " + ctesID + " ALREADY SEEN!");
             System.exit(1);
           }

           cte.append("adj:" + orient + "\n");
System.err.println("LOOKING UP IDS " + splitLine[1] + " AND " + splitLine[3]); 
System.err.println("FOUND " + idToUID.get(splitLine[1]) + " AND " + idToUID.get(splitLine[3]));
System.err.println("FINAL MAP " + eidTOiid.get(idToUID.get(splitLine[1])) + " AND " + eidTOiid.get(idToUID.get(splitLine[3])));
           if (eidTOiid.get(idToUID.get(splitLine[1])) == null || eidTOiid.get(idToUID.get(splitLine[3])) == null) {
              // skip mappings for which one of the contigs isn't in the ace file
              continue;
           }
           cte.append("nds:" + eidTOiid.get(idToUID.get(splitLine[1])) + "," + eidTOiid.get(idToUID.get(splitLine[3])) + "\n");
           cte.append("iid:" + iid + "\n");
           cte.append("sze:0\n");
           cte.append("lnk:\n");
           iid++;

           for (int i = 0; i < Integer.parseInt(splitLine[5]); i++) {
              StringBuffer ctl = new StringBuffer();
              ctl.append("{CTL\nadj:" + orient + "\nnds:" + eidTOiid.get(idToUID.get(splitLine[1])) + "," + eidTOiid.get(idToUID.get(splitLine[3])) + "\niid:" + iid + "\ntyp:O\nsze:0\nsrc:0,OVL\nstd:0\n}\n");
              cte.append(iid + "\n");
              iid++;

              System.out.print(ctl);
           }
           cte.append(".\n");
           cte.append("}\n");
           ctes.put(ctesID, cte.toString());
         } else if (line.startsWith("F")) {
         } else if (line.startsWith("I")) {
         } else if (line.startsWith("P")) {
         } else if (line.startsWith("S")) {
         } else {
           idToUID.put(splitLine[0], splitLine[1]);
System.err.println("STORING MAPPING " + splitLine[0] + " TO " + splitLine[1]);
         }
      }
      bf.close();

      if (outputCTE) {
         for (String s : ctes.keySet()) {
            System.out.print(ctes.get(s));
         }
      }
   }
}
