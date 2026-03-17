package com.example;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap; 


//Class that defines an object
class Gene{
    //attributes to store data and differentiate from other genes
    private String name;
    private String sequence;
    private List<Gene> children;
    private Gene parent;
    private double parental_distance;
    private Map<String, Map<Double, String>> mutationPattern;


    //Basic setters and getters
    public void setsequence(String sequence){
        this.sequence = sequence;
    }


    public void setMutationRate(Map<String, Map<Double, String>> mutationPattern){
        this.mutationPattern = mutationPattern;
    }

    public void setName(String name){
        this.name = name;
    }


    public String getName(){
        return name;
    }
    
     public String getSequence(){
        return sequence;
    }

     public Map<String, Map<Double, String>> getMutationPattern(){
        return mutationPattern;
    }

    public List<Gene> getChildren(){
        return children;
    }

    public Gene getParent(){
        return parent;
    }

    public double getParental_distance(){
        return parental_distance;
    }

    //Constructors
    public Gene(String name,String sequence, Map<String, Map<Double, String>> mutationPattern){
        this.name = name.replaceAll("[^A-Za-z0-9_.-]", "_");
        this.sequence = sequence;
        this.mutationPattern = mutationPattern;
        this.parental_distance = 0;
        this.parent = null;
        this.children = new ArrayList<>();
    }

    public Gene(String name,String sequence, Map<String, Map<Double, String>> mutationPattern,double parental_distance, Gene parent){
        this.name = name.replaceAll("[^A-Za-z0-9_.-]", "_");
        this.sequence = sequence;
        this.mutationPattern = mutationPattern;
        this.parental_distance = parental_distance;
        this.parent = parent;
        this.children = new ArrayList<>();
    }

    //Adds a genetic child to the children list
    public void addChild(Gene child){
        children.add(child);
    }

    //Sums all the leaves of the recursive method below into one fasta file.
    public void writeLeavesToFasta(Gene root, String path) throws IOException {

    Path outputFile = Path.of(path);
    Files.createDirectories(outputFile.getParent());
    try (BufferedWriter w = Files.newBufferedWriter(outputFile)) {
        writeLeavesRec(root, w);
    }
}

//Recursive method to get allignment fasta files text - all the leaves of the evolution
private void writeLeavesRec(Gene thisGene, BufferedWriter writer) throws IOException {

    if (thisGene == null) return;

    // leaf
    if (thisGene.getChildren().size() == 0) {

        writer.write(">");
        writer.write((thisGene.getName()));
        writer.newLine();

        String seq = thisGene.getSequence();

        for (int i = 0; i < seq.length(); i += 60) {
            int end = Math.min(i + 60, seq.length());
            writer.write(seq, i, end - i);
            writer.newLine();
        }

        return;
    }

    for (Gene child : thisGene.getChildren()) {
        writeLeavesRec(child, writer);
    }
}


//This method write this single sequence and Gene information into a fasta file
 public void toFasta(String folderPath) throws IOException {

    Path dir = Path.of(folderPath);
    Files.createDirectories(dir);

    Path fastafile = dir.resolve(name + ".fasta");

    try (BufferedWriter w = Files.newBufferedWriter(fastafile)) {

        w.write(">");
        w.write(name);
        w.newLine();

        for (int i = 0; i < sequence.length(); i += 60) {
            int end = Math.min(i + 60, sequence.length());
            w.write(sequence, i, end - i);
            w.newLine();
        }
    }
}

}