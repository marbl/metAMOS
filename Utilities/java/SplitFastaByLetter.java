import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class SplitFastaByLetter {  
   private String letter = "N";
 
   private static final String[] FASTA_ENDS = {"fasta", "fa"}; 
   private static final NumberFormat nf = new DecimalFormat("############.#");
   
   public SplitFastaByLetter() {
   }
  
   public void outputFasta(String fastaSeq, String ID) {
      if (fastaSeq.length() == 0) {
         return;
      }
      String[] split = fastaSeq.trim().toUpperCase().split(letter + "+");

      for (int i = 0; i < split.length; i++) {
         System.out.println(">" + ID + "_" + i);
         System.out.println(Utils.convertToFasta(split[i]));
      }
   }
 
   public void processFile(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FASTA_ENDS);
      
      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            outputFasta(fastaSeq.toString(), header);
            header = line.split("\\s+")[0].substring(1);
            fastaSeq = new StringBuffer();
         }
         else {
            fastaSeq.append(line);
         }
      }

      outputFasta(fastaSeq.toString(), header);
      bf.close();
   }

   public static void printUsage() {
      System.err.println("This program splits fasta records by a specified string sequence. The default sequence is N. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: SplitFastaByLetter fasta1.fasta,fasta2.fasta NNN");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      SplitFastaByLetter f = new SplitFastaByLetter();

      if (args.length > 1) {
         f.letter = args[1];
      }
      
      String[] splitLine = args[0].trim().split(",");
      for (int j = 0; j < splitLine.length; j++) {
System.err.println("Processing file " + splitLine[j]);
     	  f.processFile(splitLine[j]);
      }
   }
}
