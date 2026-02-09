package com.example;
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
        this.name = name;
        this.sequence = sequence;
        this.mutationPattern = mutationPattern;
        this.parental_distance = 0;
        this.parent = null;
        this.children = new ArrayList<>();
    }

    public Gene(String name,String sequence, Map<String, Map<Double, String>> mutationPattern,double parental_distance, Gene parent){
        this.name = name;
        this.sequence = sequence;
        this.mutationPattern = mutationPattern;
        this.parental_distance = parental_distance;
        this.parent = parent;
        this.children = new ArrayList<>();
    }


    public void addChild(Gene child){
        children.add(child);
    }

}