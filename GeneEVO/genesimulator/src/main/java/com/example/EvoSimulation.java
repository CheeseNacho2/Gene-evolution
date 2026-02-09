package com.example;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;

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

    // One instance of simulation - making one time step - it is a general method,
    // we will need to make
    // a method to calculate the mutation rates
    public Gene simulateOne(Gene SimulGene, double branch_length) {
        String evolvedGC = "";
        for (int k = 0; k < (SimulGene.getSequence()).length(); k++) {
            Random rand = new Random();
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
                evolvedGC += transitMap.get((Math.round(transition/((0.25 - 0.25*Math.exp((-4/3)*branch_length*1))*100)*10.0))/10.0); //This is to prevent errors with floating point numbers
            }

            else{
                evolvedGC += transitMap.get(transition);
            }

        }
        Gene evolvedGene = new Gene(SimulGene.getName() + "-" + Integer.toString(SimulGene.getChildren().size()), evolvedGC, SimulGene.getMutationPattern(),branch_length, SimulGene);
        SimulGene.addChild(evolvedGene);
        return evolvedGene;
    }

    // Making multiple time steps - try make it work for different time steps each
    // time
    public List<Gene> simulateMany(Gene SimuGene, double[] branch_lengths) {
        List<Gene> newGenes = new ArrayList<Gene>();
        for (int k = 0; k < branch_lengths.length; k++) {
            newGenes.add(simulateOne(SimuGene, branch_lengths[k]));
        }
        return newGenes;
    }

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
        newickform.append(gene.getName());
    }

    return newickform.toString();
}


    // Puts a gene into the database for simulation
    public Gene createRootGene(String name, String sequence, Map<String, Map<Double, String>> mutationRate) {
        rootGenes.put(name, new Gene(name, sequence, mutationRate));
        return rootGenes.get(name);
    }

    public Map<String, Map<Double, String>> JukesCantor(String sequence) {
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

    public static void main(String args[]) {
        EvoSimulation tester = new EvoSimulation();
        Gene TestingSample = new Gene("Test", "ATACGTAGCATCGATCGATCGATCGACATCGATCGATCGATCGACTGCATCGATCGACTGCATCGACTACGATCGACTAGCTAAATCGCGTCGATCGATGCTAGTCAGCTGATCGATCGATCGAGCTAGCATCGATCAGCTAG",
                tester.JukesCantor("ATACGTAGCATCGATCGATCGATCGACATCGATCGATCGATCGACTGCATCGATCGACTGCATCGACTACGATCGACTAGCTAAATCGCGTCGATCGATGCTAGTCAGCTGATCGATCGATCGAGCTAGCATCGATCAGCTAG"));
        
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
        String new_gen = tester.simulateOne(TestingSample, branch_length).getSequence();
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
        for (int i = 0; i < 5; i++){
            Random rand = new Random();
            double randomBranch = rand.nextInt(50)/10.0;
            double randomBranch2 = rand.nextInt(50)/10.0;
            Gene child_one = tester.simulateOne(evolvingGene, randomBranch);
            Gene child_two = tester.simulateOne(evolvingGene, randomBranch2);
            
            if (i % 2 == 0){
                evolvingGene = child_one;
            }
            else{
                evolvingGene = child_two;
            }
        }

        tester.newick(TestingSample.getChildren().get(0));
        System.out.println(tester.getNewickFormat().get(TestingSample.getChildren().get(0)));

    }
}
