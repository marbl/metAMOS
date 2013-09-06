import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class convertFastaAndQualToFastq {  
   private static final NumberFormat nf = new DecimalFormat("############.#");
   private static char OFFSET_CHAR = '!'; 
   private static final String[] FILE_SUFFIX = {"seq", "fasta", "fna", "fa"};

   public convertFastaAndQualToFastq() {
   }

   public void processFasta(String inputFile, String quals) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FILE_SUFFIX);
      if (bf == null) { return; }
      BufferedReader qual = (quals == null || "null".equalsIgnoreCase(quals) ? null : Utils.getFile(quals, "qual"));

      String line = null;
      StringBuilder fastaSeq = new StringBuilder();
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            if (header.length() > 0) {
               outputFasta(fastaSeq.toString(), Utils.encodeQualRecord(qual, fastaSeq.length(), OFFSET_CHAR), header, "@", "+", false);
            }
            header = line.substring(1);
            fastaSeq = new StringBuilder();
         }
         else {
            fastaSeq.append(line);
         }
      }

      if (fastaSeq.length() != 0) { 
               outputFasta(fastaSeq.toString(), Utils.encodeQualRecord(qual, fastaSeq.length(), OFFSET_CHAR), header, "@", "+", false);
      }
      if (qual != null) qual.close();
      bf.close();
   }

   public void outputFasta(String fastaSeq, String qualSeq, String ID, String fastaSeparator, String qualSeparator, boolean convert) {
      if (fastaSeq.length() == 0) {
         return;
      }

      if (qualSeq != null && qualSeq.length() != fastaSeq.length()) {
         System.err.println("Error length of sequences and fasta for id " + ID + " aren't equal fasta: " + fastaSeq.length() + " qual: " + qualSeq.length());
         System.exit(1);
      }

         System.out.println(fastaSeparator + ID);
         System.out.println((convert == true ? Utils.convertToFasta(fastaSeq) : fastaSeq));

         if (qualSeq != null) {
            System.out.println(qualSeparator + ID);
            System.out.println((convert == true ? Utils.convertToFasta(qualSeq) : qualSeq));
         }
      }

   public static void printUsage() {
      System.err.println("This program converts and fasta and qual file to fastq.");
      System.err.println("Example usage: convertFastaAndQualToFastq fasta.fasta fasta.qual encoding");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}
      if (args.length > 2) { convertFastaAndQualToFastq.OFFSET_CHAR = args[2].charAt(0); }

      convertFastaAndQualToFastq f = new convertFastaAndQualToFastq();

      f.processFasta(args[0], (args.length >= 2 ? args[1]: null));
   }
}
