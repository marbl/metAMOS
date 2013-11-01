import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;
import java.io.PrintStream;
import java.io.File;

public class convertFastqToFasta {  
  private static final NumberFormat nf = new DecimalFormat("############.#");
  private static char OFFSET_CHAR = '!';
  public static final Pattern splitBySpaces = Pattern.compile("\\s+");
  public static final String[] FILE_SUFFIX = {"seq", "fastq", "fq"};
 
   public convertFastqToFasta() {
   }

   public void processFasta(String inputFile, String outputFasta, String outputQual) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FILE_SUFFIX);
      if (bf == null) { return; }
      PrintStream fastaOut = new PrintStream(new File(outputFasta));
      PrintStream qualOut = new PrintStream(new File(outputQual));

      String line = null;
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         fastaOut.println(">"+line.substring(1, line.length()));
         line = bf.readLine();
         fastaOut.println(Utils.convertToFasta(line));
         line = bf.readLine();
         qualOut.println(">"+line.substring(1, line.length()));
         line = bf.readLine();
         qualOut.println(Utils.convertToFasta(Utils.decodeQualRecord(line, line.length(), OFFSET_CHAR)));
      }
      bf.close();
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: convertFastqToFasta in.fq out.fasta out.qual encoding");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 3) { printUsage(); System.exit(1);}
      if (args.length > 3) { convertFastqToFasta.OFFSET_CHAR = args[3].charAt(0); }

      convertFastqToFasta f = new convertFastqToFasta();
      f.processFasta(args[0], args[1], args[2]);
   }
}
