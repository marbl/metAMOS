import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import java.io.File;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipEntry;

public class SizeFasta {  
  private static final NumberFormat nf = new DecimalFormat("############");
   private static final String[] suffix = {"seq", "contig", "contigs", "scaftig", "scafSeq", "fa", "fna", "fasta", "final"}; 
   private static final String[] suffixFQ = {"fastq", "fq", "txt"};

   private double totalLen = 0;
   private boolean totalOnly = false;
   private boolean ungapped = false;

   public SizeFasta() {
   }

   private int getFastaStringLength(StringBuffer fastaSeq) {
      return (ungapped == false ? fastaSeq.length() : fastaSeq.toString().replaceAll("N", "").replaceAll("n", "").replaceAll("-", "").length());
   }

   public void processFasta(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, suffix);
      
      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            if (fastaSeq.length() != 0) {
               if (totalOnly) { totalLen += getFastaStringLength(fastaSeq); }
               else { System.out.println(header + "\t" + getFastaStringLength(fastaSeq)); }
            }
            header = line.substring(1);
            fastaSeq = new StringBuffer();
         }
         else {
            fastaSeq.append(line);
         }
      }

      if (fastaSeq.length() != 0) { 
         if (totalOnly) { totalLen += getFastaStringLength(fastaSeq); }
         else { System.out.println(header + "\t" + getFastaStringLength(fastaSeq)); }
      }
      bf.close();
   }

   public void processFastq(String inputFile) throws Exception {
      BufferedReader bf = null;

      if (inputFile.endsWith("zip")) {
         ZipInputStream zi = new ZipInputStream(new FileInputStream(new File(inputFile)));
         ZipEntry ze = null;
         while ((ze = zi.getNextEntry()) != null) {
            if (ze.isDirectory() == true) {
               System.err.println("Error, only support files not directories in zip files!");
               System.exit(1);
            }
System.err.println("Zip entry is " + ze);
            bf = new BufferedReader(new InputStreamReader(zi));
            processFastq(bf);
         }
         zi.close();
      } else {
         bf = Utils.getFile(inputFile, suffixFQ);
         processFastq(bf);
         bf.close();
      }
   }

   public void processFastq(BufferedReader bf) throws Exception {
      if (bf == null) { return; }

      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";

      while ((line = bf.readLine()) != null) {
         // read four lines at a time for fasta, qual, and headers
         String ID = line.split("\\s+")[0].substring(1);
         String fasta = bf.readLine();
         String qualID = bf.readLine().split("\\s+")[0].substring(1);

         if (qualID.length() != 0 && !qualID.equals(ID)) {
            System.err.println("Error ID " + ID + " DOES not match quality ID " + qualID);
            System.exit(1);
         }
         String qualSeq = bf.readLine();
         if (totalOnly) {
            totalLen += fasta.length();
         } else { 
            System.out.println(ID + "\t" + fasta.length());
         }
      }
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: SizeFasta fasta1.fasta,fasta2.fasta");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      SizeFasta f = new SizeFasta();
      if (args.length >= 2) { f.ungapped = Boolean.parseBoolean(args[1]); }

      boolean processed = false; 
      String[] splitLine = args; //args[0].trim().split(",");
      for (int j = 0; j < splitLine.length; j++) {
         if (splitLine[j].equalsIgnoreCase("-t")) {
            f.totalOnly = true;
            continue;
         }

         processed = false;

          for (String s : SizeFasta.suffixFQ) {
             if (splitLine[j].contains(s)) {
                f.processFastq(splitLine[j]);
                processed = true;
                break;
             }
          }
          if (!processed) {
     	  for (String s : SizeFasta.suffix) {
             if (splitLine[j].contains(s)) {
                f.processFasta(splitLine[j]);
                processed = true;
                break;
             }
           }
          } 
          if (!processed){
             System.err.println("Unknown file type " + splitLine[j]);
          }
       }

       if (f.totalOnly) {
          System.out.println(nf.format(f.totalLen));
       }
   }
}
