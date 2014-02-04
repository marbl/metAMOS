import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

public class SubFile {
   public static final int BASES_IN=40;
   public static final double THRESHOLD=0.00;
   public static final boolean OUTPUT_OVERLAP = true;
   
   public static final int KEY_ID = 0;

   private class Pair {
      Pair(double id, double beg, double end) {
         this.beg = Math.min(beg,end);
         this.end = Math.max(beg,end);
         this.id = id;
      }
      Pair(double id, double beg, double end, String other) {
         this(id, beg, end);
         info = other;
      }
      public double beg;
      public double end;
      public double id;
      public String info;
   }
   
   private HashMap<String, Integer> toRead = new HashMap<String, Integer>();
   
   public void subFile(String keysFile, String inputFile, Integer idCol, Integer idCol2, boolean rev) throws Exception {
      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(keysFile)));
      String line = null;
      String prev = null;

System.err.println("BEGIN PROCESSING KEYS");

HashMap<String, ArrayList<String>> lineNum = new HashMap<String, ArrayList<String>>();
      while ((line = bf.readLine()) != null) {
         String[] splitLine = line.trim().split("\\s+");
         
         try {
               if (!toRead.containsKey(splitLine[KEY_ID])){
                  toRead.put(splitLine[KEY_ID], 0);  
               }
               //result.put(splitLine[KEY_ID] + toRead.get(splitLine[KEY_ID]).toString(), "");
               toRead.put(splitLine[KEY_ID], toRead.get(splitLine[KEY_ID])+1);               
               if (lineNum.get(splitLine[KEY_ID]) == null) {
                  lineNum.put(splitLine[KEY_ID], new ArrayList<String>());
               }
               lineNum.get(splitLine[KEY_ID]).add(splitLine[0]);
         } catch (Exception e) {
            System.err.println("Warning: Could not parse line " + line);
         }
      }
System.err.println("DONE PROCESSING KEYS");      
      bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(inputFile)));
      int count = 0;

      while ((line = bf.readLine()) != null) {
         if (line.startsWith("uid")) { System.out.println(line); continue;}
         if (line.startsWith("ID")) { continue; }
         
         if (count % 1000000 == 0) { System.err.println("PROCESSED " + count + " RECORDS FOR OUTPUT"); }
         String[] splitLine = line.trim().split("\\s+");
         int output = 0;

         if (splitLine.length <= idCol) {
        	 continue;
         }
         String id = splitLine[idCol];//.substring(0, 4);
         if (idCol2 == -1) {
            if (toRead.containsKey(id)) {
               output = toRead.get(id);               
            }
         }
         
         if (rev) {
            if (output == 0) { output++; }
            else { output = 0; }
         }

         for (int i = 0; i < output; i++) {
            System.out.println(line);
            count++;
         }         
         prev = line;
         
         count++;
      }
   }
   
   public static void main(String[] args) throws Exception {
      int idCol = 0;
      int idCol2 = -1;
      boolean rev = false;
      
      if (args.length >= 3) {
         idCol = Integer.parseInt(args[2]);
      }
      if (args.length >=4) {
         idCol2 = Integer.parseInt(args[3]);
      }
      if (args.length >= 5) {
         rev = Boolean.parseBoolean(args[4]);
      }

      SubFile f = new SubFile();
      f.subFile(args[0], args[1], idCol, idCol2, rev);
   }
}
