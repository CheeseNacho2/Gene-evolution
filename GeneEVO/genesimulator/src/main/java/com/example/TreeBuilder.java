package com.example;
import java.util.Random;


public class TreeBuilder {

    public static void binaryTree(int currentGeneration, int maxGenerations, Gene rootGene) {

        if (currentGeneration > maxGenerations) {
            return;
        }

        Gene currentGene = rootGene;

        Random rand = new Random();
        double randomBranch = rand.nextInt(50)/10.0;
        double randomBranch2 = rand.nextInt(50)/10.0;

        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch));
        binaryTree(currentGeneration + 1, maxGenerations, EvoSimulation.simulateOne(currentGene, randomBranch2));

    }
}