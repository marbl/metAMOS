import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.TreeMap;

/** 
 * Class to output default quality values given a fasta file
 * Necessary when a program expects both quality and bases but we only have baess
 */
class outputDefaultQuality {
   private static int DEFAULT_QUAL = 40;
   private static final int LINE_WIDTH = 20;

   private static String buildQuality(StringBuilder fasta) {
      StringBuilder st = new StringBuilder();

      if (fasta.length() != 0) {
             for (int i = 0; i < fasta.length(); i++) {
                if (i != 0 && i % LINE_WIDTH == 0) { st.append("\n"); }
                   if (i != 1 && ((i+1) % LINE_WIDTH == 0)) st.append(DEFAULT_QUAL); else st.append(DEFAULT_QUAL + " ");
              }
      }
      return st.toString();
   }

   public static void printUsage() {
      System.err.println("This program outputs default quality values given a fasta file. The default quality is " + DEFAULT_QUAL);
      System.err.println("Example usage: java outputDefaultQuality file.fasta <QUALITY_VALUE>");
   }

   public static void main(String[] args) throws Exception {
      if (args.length < 1) {printUsage(); System.exit(1); }
      if (args.length >= 2) { DEFAULT_QUAL = Integer.parseInt(args[1]); }

      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(args[0])));
      String line = null;
      StringBuilder totalLine = new StringBuilder();
      String header = null;

      while ((line = bf.readLine()) != null) {
        if (line.startsWith(">")) {
           if (totalLine.length() != 0) {
              System.out.println(header);
              System.out.println(buildQuality(totalLine));
           }
           header = line; 
           totalLine = new StringBuilder();
        }
        else {
           totalLine.append(line);
        }
     }
     if (totalLine.length() != 0) {
        System.out.println(header);
        System.out.println(buildQuality(totalLine));
     }
     bf.close(); 
  }
}
