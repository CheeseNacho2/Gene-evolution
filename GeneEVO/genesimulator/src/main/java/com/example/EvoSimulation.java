package com.example;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.File;   
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.apache.commons.math3.distribution.GammaDistribution;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtils;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.w3c.dom.Node;

import com.example.Gene;

import javax.swing.*;
import java.io.File;

//Class for simulation
class EvoSimulation {
    private Map<String, Gene> rootGenes;
    private Map<Gene, String> newickformat;

    // Constructor
    public EvoSimulation() {
        this.rootGenes = new HashMap<>();
        this.newickformat = new HashMap<>();
    }


    // Getters
    public Map<String, Gene> getRootGenes(){
        return rootGenes;
    }

    public Map<Gene, String> getNewickFormat(){
        return newickformat;
    }

    public void createRootGene(String name, Gene rootGene){
        rootGenes.put(name, rootGene);
    }

    public void newick(Gene gene){
        newickformat.put(gene, toNewick(gene) + ";");
    }

    //Class for getting a variable-sized array with Gamma Distributed values 
    public static double[] gammaDistributeArray(double shape, int element_count){
        double scale = 1/shape;
        GammaDistribution gamma = new GammaDistribution(shape, scale);
        double[] gamma_distributed_array = new double[element_count];
        for (int i = 0; i < element_count; i++ ){
            gamma_distributed_array[i] = gamma.sample();
        }

        return gamma_distributed_array;
    }


    //This method converts fasta files into Gene objects using the upper most sequence
    public static Gene fastafile(String filename, String model) throws IOException {

    InputStream filestream = EvoSimulation.class.getClassLoader().getResourceAsStream(filename);

        if (filestream == null) {
            throw new RuntimeException("File not found in resources!");
    }
    if (!"Jukes-Cantor".equals(model) && !"HKY".equals(model)) {
        model = "Jukes-Cantor";
        System.out.println("Gene was created using Jukes-Cantor model due to unknown model input.");
    }
     try (BufferedReader reader = new BufferedReader(new InputStreamReader(filestream))) {

        StringBuilder sequence = new StringBuilder();
        String name = null;
        String readingline;

        while ((readingline = reader.readLine()) != null) {

            if (readingline.startsWith(">")) {

                if (name == null) {

                    int dot = readingline.indexOf('.');

                    if (dot == -1 || dot + 2 >= readingline.length()) {
                        throw new IOException("Wrong header format");
                    }

                    name = readingline.substring(dot + 3, readingline.indexOf(","));

                } else {
                    break;
                }

            } else {
                sequence.append(readingline.trim());
            }
        }

        if (name == null) {
            throw new IOException("Wrong file format.");
        }

        String seq = sequence.toString();

        if ("Jukes-Cantor".equals(model)) {
            return new Gene(name, seq, JukesCantor(seq));
        } else {
            return new Gene(name, seq, HKY85(seq));
        }
    }
}
    // One instance of simulation - making one time step - it is a general method,
    // we will need to make
    // a method to calculate the mutation rates
    public static Gene simulateOne(Gene SimulGene, double branch_length) {
        StringBuilder evolvedGC = new StringBuilder();
        double og_branch_length = branch_length;
        double[] gamma_distributed_array = gammaDistributeArray(0.6, SimulGene.getSequence().length());
        Random rand = new Random();
        for (int k = 0; k < (SimulGene.getSequence()).length(); k++) {
            branch_length = gamma_distributed_array[k];
            int randomValue = rand.nextInt(101);
            double transition = randomValue;
            Map<Double, String> transitMap = ((SimulGene.getMutationPattern())
                    .get(Character.toString((SimulGene.getSequence()).charAt(k))));
            Double prevKey = 101.0;
            for (Double probNumber : transitMap.keySet()) {
                if (probNumber < 100) {
                    probNumber *= (0.25 - 0.25*Math.exp((-4/3)*branch_length*1)) * 100.0;
                }
                if (randomValue < probNumber && prevKey > probNumber) {
                    transition = probNumber;
                    prevKey = probNumber;
                }
            }

            //add error handling for null values
            if(transition != 100.0){
                evolvedGC.append(transitMap.get((Math.round(transition/((0.25 - 0.25*Math.exp((-4/3)*branch_length*1))*100)*10.0))/10.0)); //This is to prevent errors with floating point numbers
            }

            else{
                evolvedGC.append(transitMap.get(transition));
            }

        }
        Gene evolvedGene = new Gene(SimulGene.getName() + "-" + (char)('A' + (SimulGene.getChildren().size())), evolvedGC.toString(), SimulGene.getMutationPattern(),og_branch_length, SimulGene);
        SimulGene.addChild(evolvedGene);
        return evolvedGene;
    }

    // Making multiple time steps - try make it work for different time steps each
    // time
    public List<Gene> simulateMany(Gene SimuGene, double[] branch_lengths, String model) {
        List<Gene> newGenes = new ArrayList<Gene>();
        for (int k = 0; k < branch_lengths.length; k++) {
            newGenes.add(simulateOne(SimuGene, branch_lengths[k]));
        }
        return newGenes;
    }

    //Method to create a string that represents a tree in Newick format - without the root
    public static String toNewick(Gene gene) {
    if (gene.getChildren().isEmpty()) {
        return gene.getName() + ":" + String.format(Locale.UK, "%.3f", gene.getParental_distance());
    }

    StringBuilder newickform = new StringBuilder();
    newickform.append("(");

    for (int i = 0; i < gene.getChildren().size(); i++) {
        newickform.append(toNewick(gene.getChildren().get(i)));

        if (i < gene.getChildren().size() - 1) {
            newickform.append(",");
        }
    }

    newickform.append(")");

    if (gene.getName() != null && !gene.getName().isEmpty()) {
        newickform.append(//gene.getName() + 
        ":" + String.format(Locale.UK, "%.3f", gene.getParental_distance()));
    }

    return newickform.toString();
}


    // Puts a gene into the database for simulation
    public Gene createRootGene(String name, String sequence, Map<String, Map<Double, String>> mutationRate) {
        rootGenes.put(name, new Gene(name, sequence, mutationRate));
        return rootGenes.get(name);
    }

    public static Map<String, Map<Double, String>> JukesCantor(String sequence) {
        Map<String, Map<Double, String>> mutationRates = new HashMap<>();
        double JC_rate = 1.0;
        String[] nucleotides = { "A", "C", "T", "G" };
        for (String x : nucleotides) {
            mutationRates.put(x, new HashMap<>());
            int probability_handler = 1;
            for (String y : nucleotides) {
                if (x.equals(y)) {
                    mutationRates.get(x).put(1.0 * 100, y);
                } else {
                    mutationRates.get(x).put(probability_handler * JC_rate * 1.0, y);
                    probability_handler++;
                }
            }
        }

        return mutationRates;
    }

    //HKY method for the initial rate matrix - it is purely a rate matrix and is not used in actual simulation due to high complexity of matrix exponentials
    public static Map<String, Map<Double, String>> HKY85(String sequence){
        Map<String, Map<Double, String>> mutationRates = new HashMap<>();
        Map<String, Double> base_frequencies = new HashMap<>();
        double transition_bias = 5.0;
         String[] nucleotides = { "A", "C", "T", "G" };
        base_frequencies.put("A", 0.0);
        base_frequencies.put("C", 0.0);
        base_frequencies.put("T", 0.0);
        base_frequencies.put("G", 0.0);
        for(int x = 0; x < sequence.length() - 1; x++){
            base_frequencies.put(String.valueOf(sequence.charAt(x)), base_frequencies.get(String.valueOf(sequence.charAt(x))) + 1);
        }
        base_frequencies.put("A",base_frequencies.get("A")/sequence.length());
        base_frequencies.put("C",base_frequencies.get("C")/sequence.length());
        base_frequencies.put("T",base_frequencies.get("T")/sequence.length());
        base_frequencies.put("G",base_frequencies.get("G")/sequence.length());
         for (String x : nucleotides) {
            mutationRates.put(x, new HashMap<>());
            double probability_handler = 100.0;
            
            for (String y : nucleotides) {
                if (x.equals(y)) {
                    mutationRates.get(x).put(1.0 * 100, y);
                } 
                else if(((x == ("A") || x == "G" ) && (y == "A" || y == "G")) || ((x == ("C") || x == "T" ) && (y == "C" || y == "T"))){
                    probability_handler += base_frequencies.get(y) * transition_bias;
                    mutationRates.get(x).put(probability_handler * 100, y);
                }
                else {
                    probability_handler += base_frequencies.get(y);
                    mutationRates.get(x).put(probability_handler * 100, y);
                    
                }
                
            }
             mutationRates.get(x).put(probability_handler, x);
            
        }
        
        return mutationRates;
    }


//Method for simulating via HKY85 method - uses equations for each transition
public static Gene simulateOne(Gene evolving_gene, double branch_length, double transition_bias){
    StringBuilder evolvedGC = new StringBuilder();
    Map<String, Map<Double, String>> mutationRates = new HashMap<>();
    Map<String, Double> base_frequencies = new HashMap<>();
    String[] nucleotides = {"A", "C", "T", "G"};
    for (String nucleotide : nucleotides){
    base_frequencies.put(nucleotide, 0.0);
    }
    for(int x = 0; x < evolving_gene.getSequence().length() - 1; x++){
        base_frequencies.put(String.valueOf(evolving_gene.getSequence().charAt(x)), base_frequencies.get(String.valueOf(evolving_gene.getSequence().charAt(x))) + 1);
    }
    for (String nucleotide : nucleotides){
     base_frequencies.put(nucleotide,base_frequencies.get(nucleotide)/evolving_gene.getSequence().length());
    }

    double beta = 1/(2*(base_frequencies.get("A") + base_frequencies.get("G")) * (base_frequencies.get("C") + base_frequencies.get("T")) + 
                  2*transition_bias*((base_frequencies.get("A")*base_frequencies.get("G")) +(base_frequencies.get("C")*base_frequencies.get("T"))));

    for (String x : nucleotides) {
            mutationRates.put(x, new HashMap<>());
            double probability_handler = 0.0;
            
            for (String y : nucleotides) {
                String[] transversion_nucleotides = new String[2];
                if (x.equals(y)) {
                    mutationRates.get(x).put(1.0 * 100, y);
                    System.out.println("unchanged");
                } 
                else if(((x == ("A") || x == "G" ) && (y == "A" || y == "G")) || ((x == ("C") || x == "T" ) && (y == "C" || y == "T"))){
                    if(x == ("A") || x == "G" ){
                        transversion_nucleotides[0] = "C";
                        transversion_nucleotides[1] = "T";
                    }
                    else
                    {
                        transversion_nucleotides[0] = "A";
                        transversion_nucleotides[1] = "G";
                    }
                    double transition_rate = (base_frequencies.get(y)*
                    (base_frequencies.get(x) + base_frequencies.get(y) + 
                    (base_frequencies.get(transversion_nucleotides[0])+base_frequencies.get(transversion_nucleotides[1]))*Math.exp(-beta*branch_length))
                    - base_frequencies.get(y)*Math.exp(-(1+(base_frequencies.get(x) + base_frequencies.get(y))*(transition_bias - 1))*beta*branch_length))/
                    (base_frequencies.get(x) + base_frequencies.get(y));
                    mutationRates.get(x).put(transition_rate * 100 + probability_handler, y);
                    probability_handler += transition_rate * 100;
                    System.out.println(transition_rate + "  transition from " + x + " to " + y);
                }
                else {
                    double transition_rate = base_frequencies.get(y)*(1.0 - Math.exp(-beta*branch_length));
                    mutationRates.get(x).put(transition_rate * 100 + probability_handler, y);
                    probability_handler += transition_rate * 100;
                    System.out.println(transition_rate + "  transversion from"+ x + " to " + y);
                }
                
            }
        }
            Random rand = new Random();
            for (int k = 0; k < (evolving_gene.getSequence()).length(); k++) {
            
            int randomValue = rand.nextInt(101);
            double transition = randomValue;
            Map<Double, String> transitMap = (mutationRates.get(Character.toString((evolving_gene.getSequence()).charAt(k))));
            Double prevKey = 101.0;
            for (Double probNumber : transitMap.keySet()) {
                if (randomValue < probNumber && prevKey > probNumber) {
                    transition = probNumber;
                    prevKey = probNumber;
                }
            }
            evolvedGC.append(transitMap.get(transition));
        } 
        
    
    
    Gene evolvedGene = new Gene(evolving_gene.getName() + "-" + (char)('A' + (evolving_gene.getChildren().size())),
    evolvedGC.toString(), evolving_gene.getMutationPattern(),branch_length, evolving_gene);
    //System.out.println(evolvedGC);
    evolving_gene.addChild(evolvedGene);
    return evolvedGene;
    }

    //Copied JFReeChart library for charting
    public static void createXYChart(String chartTitle, double[] x_coordinate, double[] y_coordinate, String YLabel, String XLabel) {
        XYSeries series = new XYSeries(YLabel);

        for (int i = 0; i < x_coordinate.length; i++) {
            series.add(x_coordinate[i], y_coordinate[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                chartTitle,
                XLabel,
                YLabel,
                dataset
        );

        // Chart displaying
        JFrame frame = new JFrame(chartTitle);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String args[]) throws IOException{
        EvoSimulation tester = new EvoSimulation();
        Gene TestingSample = new Gene("Test", "ATACGTAGCATCGATCGATCGATCGACATCGATCGATCGATCGACTGCATCGATCGACTGCATCGACTACGATCGACTAGCTAAATCGCGTCGATCGATGCTAGTCAGCTGATCGATCGATCGAGCTAGCATCGATCAGCTAG",
                JukesCantor("ATACGTAGCATCGATCGATCGATCGACATCGATCGATCGATCGACTGCATCGATCGACTGCATCGACTACGATCGACTAGCTAAATCGCGTCGATCGATGCTAGTCAGCTGATCGATCGATCGAGCTAGCATCGATCAGCTAG"));
        
        //MANUAL BRANCH INPUT
        //Scanner in_scanner = new Scanner(System.in);
        //System.out.println("Enter a branch length");
        //double branch_length = in_scanner.nextDouble();

        tester.createRootGene(TestingSample.getName(), TestingSample);

        //
        /*double[] b_lengths = new double[1001];
        double[] max_likelihoods = new double[1001];
        double[] hamming_distances = new double[1001];

        for (double k = 0; k < 10; k += 0.01){
        double branch_length = k;
        b_lengths[(int) Math.round(k*100)] = Math.round(k*100)/100.0;
        String new_gen = simulateOne(TestingSample, branch_length).getSequence();
        String root_gen = TestingSample.getSequence();
        double hamming_distance = 0;
            for (int i = 0; i < (root_gen.length()); i++){
                //System.out.println(i);
                if (new_gen.charAt(i) != (root_gen.charAt(i))){
                    hamming_distance += 1.0;
                }
            }
        max_likelihoods[(int) Math.round(k*100)] = -0.75 * Math.log(1 - ((4.0/3.0) * (hamming_distance/root_gen.length())));
        hamming_distances[(int) Math.round(k*100)] = hamming_distance/root_gen.length();
        //System.out.println(new_gen);
        //System.out.println(root_gen);
        }
        createXYChart("Maximum Likelihood graph", b_lengths, max_likelihoods, "Maximum Likelihood", "Branch length");
        createXYChart("Hamming distance", b_lengths, hamming_distances, "Hamming distance", "Branch Length");
        */


        //in_scanner.close();

        Gene evolvingGene = TestingSample;
        Random rand = new Random();
        /*for (int i = 0; i < 5; i++){
            double randomBranch = rand.nextInt(50)/10.0;
            double randomBranch2 = rand.nextInt(50)/10.0;
            Gene child_one = simulateOne(evolvingGene, randomBranch);
            Gene child_two = simulateOne(evolvingGene, randomBranch2);
            
            if (i % 2 == 0){
                evolvingGene = child_one;
            }
            else{
                evolvingGene = child_two;
            }
        }
        */
        
        /*TreeBuilder.binaryTree( 0, 5, TestingSample);
        tester.newick(TestingSample);
        System.out.println(tester.getNewickFormat().get(TestingSample));
        */
        simulateOne(TestingSample, 3, 5);
        Gene influenza = fastafile("sequences.fasta", "Jukes-Cantor");
        System.out.println(influenza.getName());
        System.out.println(influenza.getSequence());
        Gene Sars_Cov_19 = fastafile("sequence (1).fasta", "HKY");
        System.out.println(Sars_Cov_19.getName());
        System.out.println(Sars_Cov_19.getSequence());
        
        Gene streptococcus = fastafile("sequence (2).fasta", "HKY");
        System.out.println(streptococcus.getName());
        //System.out.println(streptococcus.getSequence());

        tester.createRootGene(influenza.getName(), influenza);
        tester.createRootGene(Sars_Cov_19.getName(), Sars_Cov_19);
        tester.createRootGene(streptococcus.getName(), streptococcus);

        TreeBuilder.binaryTree( 0, 6, influenza);
        tester.newick(influenza);
        System.out.println(tester.getNewickFormat().get(influenza));

        TreeBuilder.binaryTree( 0, 6, Sars_Cov_19);
        tester.newick(Sars_Cov_19);
        System.out.println(tester.getNewickFormat().get(Sars_Cov_19));

        TreeBuilder.binaryTree( 0, 6, streptococcus);
        tester.newick(streptococcus);
        System.out.println(tester.getNewickFormat().get(streptococcus));
        
    }
}
