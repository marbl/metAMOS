import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.io.File;

public class GetFastaStats {  
   private static final int MIN_GAP_SIZE 			= 20;
   private static int MIN_LENGTH   			= 2000;
   private static final int CONTIG_AT_INITIAL_STEP 	= 1000000;
   
   private static final NumberFormat nf = new DecimalFormat("############.#");
   private static DecimalFormat largeInts = new DecimalFormat();
   private static DecimalFormatSymbols dfs = new DecimalFormatSymbols();
   private static final int MAX_DEFAULT_STRING = 4000000;
   
   private class ContigAt {
	   ContigAt(long currentBases) {
		   this.count = this.len = 0;
		   this.totalBP = 0;
		   this.goal = currentBases;
	   }
	   
	   public int count;
	   public long totalBP;
	   public int len;
	   public long goal;
   }

   String splitLetter = null;   
   boolean baylorFormat = false;
   boolean oldStyle = false;
   boolean n50only = false;
   boolean storeCtgs = false;
   HashMap<String, Integer> fastaLens = new HashMap<String, Integer>();
   HashMap<String, String> ctgs = new HashMap<String, String>();
   int totalCount = 0;
   int genomeSize = 0;
   
   public GetFastaStats(boolean useBaylor, boolean old, boolean n50, int size, String split) {
	   baylorFormat = useBaylor;
	   oldStyle = old;
           n50only = n50;
           genomeSize = size;
           splitLetter = split;

           dfs.setGroupingSeparator(',');
           largeInts.setDecimalFormatSymbols(dfs);
   }

   public void processLens(String inputFile) throws Exception {
      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(inputFile)));

      String line = null;
      StringBuffer counterStr = new StringBuffer(10);
      int counter = 1;

      while ((line = bf.readLine()) != null) {
         String[] splitLine = line.trim().split("\\s+");

         try {
         counterStr.append(counter);
         fastaLens.put(counterStr.toString(), Integer.parseInt(splitLine[0]));
         counter++;
         counterStr.delete(0, counterStr.length());
         } catch (NumberFormatException e) {
            System.err.println("Number format exception " + e.getMessage() + " while parsing line " + line);
         }
      }
      bf.close();
   }
   
   public void processFile(String inputFile) throws Exception {
      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(inputFile)));
      
      String line = null;
      StringBuffer fastaSeq = new StringBuffer(MAX_DEFAULT_STRING);
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
        	 String fastaString = fastaSeq.toString().replaceAll("-", "");
        	 fastaString = fastaString.toString().replaceAll("\\.", "");

        	 //String[] split = fastaString.trim().split("N+");
        	 //String[] split = { fastaString.replaceAll("N", "").trim() };
String[] split = { "null" };
if (splitLetter == null) {
        	 split[0] =  fastaString.trim(); 
} else {
      split = fastaString.trim().split(splitLetter + "+");
}
//System.err.println("SPLIT one ctg of length " + fastaString.length() + " INTO " + split.length);

        	 for (int i = 0; i < split.length; i++) {
        		 if (split[i].length() != 0) {
        			 fastaLens.put(header + "_" + i,split[i].length());               
        			 if (storeCtgs)
        				 ctgs.put(header + "_" + i, split[i].toString());
        		 }
        	}
            //header = line.replace(">contig_", "").trim();
            header = line.trim().split("\\s+")[0].replace(">", "").replaceAll("ctg", "").trim();
            header = line.trim().split("\\s+")[0].replace(">", "").replaceAll("scf", "").trim();
            header = line.trim().split("\\s+")[0].replace(">", "").replaceAll("deg", "").trim();
            fastaSeq.setLength(0); // = new StringBuffer();
         }
         else {
            fastaSeq.append(line);
         }
      }

      String fastaString = fastaSeq.toString().replaceAll("-", "");
      fastaString = fastaString.toString().replaceAll("\\.", "");

      String[] split = { "null" };
      if (splitLetter == null) {
                 split[0] = fastaString.trim();
      } else {
         split = fastaString.trim().split(splitLetter + "+");
      }

      for (int i = 0; i < split.length; i++) {
    	  if (fastaSeq.length() != 0) {
    		  fastaLens.put(header + "_" + i,split[i].length());
         
    		  if (storeCtgs)
    			  ctgs.put(header + "_" + i, split[i].toString());
    	  }
      }
      bf.close();
   }
   
   public String toString(boolean outputHeader, String title) {
      StringBuffer st = new StringBuffer();
      int max = Integer.MIN_VALUE;
      int min = Integer.MAX_VALUE;
      long total = 0;
      int count = 0;
      
      int n10 = 0;
      int n25 = 0;
      int n50 = 0;
      int n75 = 0;
      int n95 = 0;
      
      int n10count = 0;
      int n25count = 0;
      int n50count = 0;
      int n75count = 0;
      int n95count = 0;

      int totalOverLength = 0;
      long totalBPOverLength = 0;
      double eSize = 0;

      if (fastaLens.size() == 0) {
         System.err.println("No contigs");
         return "";
      }

      for (String s : fastaLens.keySet()) {
         int len = fastaLens.get(s);
        
         if (oldStyle == true && len <= MIN_LENGTH) {
        	 continue;
         }
         
         if (len > max) { max = len; }
         if (len < min) { min = len; }
         total += len;
         count++;
         
         if (len > MIN_LENGTH) {
        	 totalOverLength++;
        	 totalBPOverLength += len;
         }
         eSize += Math.pow(len, 2);
      }
      eSize /= (genomeSize == 0 ? totalBPOverLength : genomeSize); 
      
      // get the goal contig at X bases (1MBp, 2MBp)
      ArrayList<ContigAt> contigAtArray = new ArrayList<ContigAt>();
      if (baylorFormat == true) {
    	  contigAtArray.add(new ContigAt( 1 * CONTIG_AT_INITIAL_STEP));
    	  contigAtArray.add(new ContigAt( 2 * CONTIG_AT_INITIAL_STEP));
    	  contigAtArray.add(new ContigAt( 4 * CONTIG_AT_INITIAL_STEP));
    	  contigAtArray.add(new ContigAt(10 * CONTIG_AT_INITIAL_STEP));
      } 
      else {      
	      long step = CONTIG_AT_INITIAL_STEP;
	      long currentBases = 0;
	      while (currentBases <= total) {
	    	  if ((currentBases / step) >= 10) {
	    		  step *= 10;
	    	  }
	    	  currentBases += step;
	    	  contigAtArray.add(new ContigAt(currentBases));
	      }
      }
      ContigAt[] contigAtVals = contigAtArray.toArray(new ContigAt[0]);

      Integer[] vals = fastaLens.values().toArray(new Integer[0]);
      Arrays.sort(vals);

      long sum = 0;
      double median = 0;
      int medianCount = 0;
      int numberContigsSeen = 0;
      int currentValPoint = 0;
      for (int i = vals.length - 1; i >= 0; i--) {
    	  if (((int) (count / 2)) == i) {
    		  median += vals[i];
    		  medianCount++;
    	  }
    	  else if (count % 2 == 0 && ((((int) (count / 2)) + 1) == i)) {
    		  median += vals[i];
    		  medianCount++;
    	  }
    	  
         sum += vals[i];

         // calculate the bases at
         if (currentValPoint < contigAtVals.length && sum >= contigAtVals[currentValPoint].goal && contigAtVals[currentValPoint].count == 0) {
        	 contigAtVals[currentValPoint].count = numberContigsSeen;
        	 contigAtVals[currentValPoint].len = vals[i];
        	 contigAtVals[currentValPoint].totalBP = sum;
        	 currentValPoint++;
         }
         // calculate the NXs
         if (sum / (double)(genomeSize == 0 ? total : genomeSize) >= 0.1 && n10count == 0) {
        	 n10 = vals[i];
        	 n10count = vals.length - i;
        	 
         }
         if (sum / (double)(genomeSize == 0 ? total : genomeSize) >= 0.25 && n25count == 0) {
        	 n25 = vals[i];
        	 n25count = vals.length - i;
        	 
         }
         if (sum / (double)(genomeSize == 0 ? total : genomeSize) >= 0.5 && n50count == 0) {
        	 n50 = vals[i];
        	 n50count = vals.length - i;
        	 
         }
         if (sum / (double)(genomeSize == 0 ? total : genomeSize) >= 0.75 && n75count == 0) {
                 n75 = vals[i];
                 n75count = vals.length - i;
         }
         if (sum / (double)(genomeSize == 0 ? total : genomeSize) >= 0.95 && n95count == 0) {
        	 n95 = vals[i];
        	 n95count = vals.length - i;
        	 
         }
         
         numberContigsSeen++;
      }
      if (medianCount != 1 && medianCount != 2) {
    	  System.err.println("ERROR INVALID MEDIAN COUNT " + medianCount);
    	  System.exit(1);
      }
      
      if (oldStyle == true) {
          if (n50only == true) {
             st.append(title + "\t" + largeInts.format(n50));
          } else {
             st.append("GenomeSize: " + genomeSize + "\n");
             st.append("NumReads: " + totalCount + "\n");
             st.append("Total units: " + count + "\n");
             st.append("BasesInFasta: " + total + "\n");
             st.append("Min: " + largeInts.format(min) + "\n");
             st.append("Max: " + largeInts.format(max) + "\n");
             st.append("N10: " + largeInts.format(n10) + " COUNT: " + n10count + "\n");
             st.append("N25: " + largeInts.format(n25) + " COUNT: " + n25count + "\n");
             st.append("N50: " + largeInts.format(n50) + " COUNT: " + n50count + "\n");
             st.append("N75: " + largeInts.format(n75) + " COUNT: " + n75count + "\n");
             st.append("N95: " + largeInts.format(n95) + " COUNT: " + n95count + "\n");    	  
             st.append("E-SIZE: " + nf.format(eSize) + "\n");
          }
      } else {
	      if (outputHeader) {
	    	  st.append("Assembly");
	    	  st.append(",Unit Number");
		      st.append(",Unit Total BP");
		      st.append(",Number Units > " + MIN_LENGTH);
		      st.append(",Total BP in Units > " + MIN_LENGTH);
		      st.append(",Min");
		      st.append(",Max");
		      st.append(",Average");
		      st.append(",Median");
		      
		      for (int i = 0; i < contigAtVals.length; i++) {
		    	  if (contigAtVals[i].count != 0) {
		    		  st.append(",Unit At " + nf.format(contigAtVals[i].goal) + " Unit Count," + /*"Total Length," + */ " Actual Unit Length" );
		    	  }
		      }
		      st.append("\n");
	      }
	      
	      st.append(title);
	      st.append("," + nf.format(count));
	      st.append("," + nf.format(total));
	      st.append("," + nf.format(totalOverLength));
	      st.append("," + nf.format(totalBPOverLength));
	      st.append("," + nf.format(min));
	      st.append("," + nf.format(max));
	      st.append("," + nf.format((double)total / count));
	      st.append("," + nf.format((double)median / medianCount));
	      
	      for (int i = 0; i < contigAtVals.length; i++) {
	    	  if (contigAtVals[i].count != 0) {
	    		  st.append("," + contigAtVals[i].count + "," + /*contigAtVals[i].totalBP + "," +*/ contigAtVals[i].len);
	    	  }
	      }
      }
      
      return st.toString();
   }
   
   public static void printUsage() {
      System.err.println("This program computes total bases in fasta sequences, N50, min, and max.");

   }
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      /*
      if (args.length > 1) {
         f.storeCtgs = true;
      }
      */           
      
      boolean useBaylorFormat = false;
      boolean oldStyle = false;
      boolean n50only = false;
      int genomeSize = 0;
      int initialVal = 0;
      String splitByLetter = null;
      for (int i = 0; i < args.length; i++) {
         if (args[i].startsWith("-")) {
    	  if (args[i].trim().equalsIgnoreCase("-b")) {
    		  useBaylorFormat = true;
    	  }
    	  else if (args[i].trim().equalsIgnoreCase("-o")) {
    		  oldStyle = true;
    	  }
          else if (args[i].trim().equalsIgnoreCase("-n50")) {
                  n50only = true;
                  oldStyle = true;
          }
    	  else if (args[i].trim().equalsIgnoreCase("-genomeSize")) {
             String param = args[++i];
             try {
                genomeSize = Integer.parseInt(param);
             } catch (NumberFormatException e) {
                // assume it's a fasta file to size
                SizeFasta s = new SizeFasta();
                genomeSize = s.sizeSingleRecord(param);
             } finally {
                initialVal++;
             }
          } else if (args[i].trim().equalsIgnoreCase("-min")) {
             GetFastaStats.MIN_LENGTH = Integer.parseInt(args[++i]);
              initialVal++;
System.err.println("Set min to be " + GetFastaStats.MIN_LENGTH);
          } else if (args[i].trim().equalsIgnoreCase("-split")) {
             splitByLetter = args[++i];
             initialVal++;
          } else {
    		  System.err.println("Unknown parameter "+ args[0] + " specified, please specify -b for Baylor-style output.");
    		  System.exit(1);
    	  }
          initialVal++;
         }
      }
      for (int i = initialVal; i < args.length; i++) {
    	  String assemblyTitle = new File(args[i].trim()).getName();
          GetFastaStats f = new GetFastaStats(useBaylorFormat, oldStyle, n50only, genomeSize, splitByLetter);
          String[] splitLine = args[i].trim().split(",");

//System.err.println("The input file is " + args[i] + " AND SPLIT INTO " + splitLine.length);
          for (int j = 0; j < splitLine.length; j++) {
System.err.println("Processing file " + splitLine[j]);
                  if (splitLine[j].endsWith("lens")) { 
                     f.processLens(splitLine[j]); 
                  } else {
        	     f.processFile(splitLine[j]);
                  }
    	  }
          System.out.println(f.toString(i == initialVal, assemblyTitle));
      }            
   }
}
