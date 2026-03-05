package com.example;
import java.io.IOException;
import java.util.Random;


public class TreeBuilder {

    public static void binaryTree(int currentGeneration, int maxGenerations, Gene rootGene) throws IOException {

        if (currentGeneration > maxGenerations) {
            return;
        }

        Gene currentGene = rootGene;

        Random rand = new Random();
        double randomBranch = rand.nextDouble(1.5);
        double randomBranch2 = rand.nextDouble(1.5);

        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch));
        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch2));

    }

    public static void binaryTree(int currentGeneration, int maxGenerations, Gene rootGene, double bias) throws IOException {

        if (currentGeneration > maxGenerations) {
            return;
        }

        Gene currentGene = rootGene;

        Random rand = new Random();
        double randomBranch = rand.nextDouble(1.5);
        double randomBranch2 = rand.nextDouble(1.5);

        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch,bias));
        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch2,bias));

    }
}