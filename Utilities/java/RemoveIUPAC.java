import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class RemoveIUPAC {  
   private static final String[] FASTA_ENDS = {"contig", "final", "fasta", "fa"}; 
   private static final NumberFormat nf = new DecimalFormat("############.#");
   
   public RemoveIUPAC() {
   }
  
   public void outputFasta(String fastaSeq, String ID) {
      if (fastaSeq.length() == 0) {
         return;
      }
      fastaSeq = fastaSeq.toUpperCase();
      StringBuilder st = new StringBuilder();
      for (int i = 0; i < fastaSeq.length(); i++) {
         if (fastaSeq.charAt(i) != 'A' && fastaSeq.charAt(i) != 'C' && fastaSeq.charAt(i) != 'G' && fastaSeq.charAt(i) != 'T' && fastaSeq.charAt(i) != 'N') { 
            st.append('N');
         } else { 
            st.append(fastaSeq.charAt(i));
         }
      }
      System.out.println(">" + ID);
      System.out.println(Utils.convertToFasta(st.toString()));
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
      System.err.println("Example usage: RemoveIUPAC fasta1.fasta,fasta2.fasta NNN");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      RemoveIUPAC f = new RemoveIUPAC();

      String[] splitLine = args[0].trim().split(",");
      for (int j = 0; j < splitLine.length; j++) {
System.err.println("Processing file " + splitLine[j]);
     	  f.processFile(splitLine[j]);
      }
   }
}
