package edu.temple.gymminder.models;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Created by rober_000 on 1/31/2017.
 */

public class Exercise implements Serializable {
    public String name;
    public int sets;
    public int reps;
    public ArrayList<Integer> completed = null;
    public ArrayList<Float> accelerationStream = null;
    public int setsDone;

    public Exercise() {

    }

    public Exercise(String name, int sets, int reps) {
        this.name = name;
        this.sets = sets;
        this.reps = reps;
    }

    public Exercise(String name, int sets, int reps, ArrayList<Integer> completed, int setsDone) {
        this.name = name;
        this.sets = sets;
        this.reps = reps;
        this.completed = completed;
        this.setsDone = setsDone;
    }

    public Exercise(Exercise exercise) {
        name = exercise.name;
        sets = exercise.sets;
        reps = exercise.reps;
        initActive();
    }

    public void initActive() {
        completed = new ArrayList<>(sets);
        for (int i = 0; i < sets; i++) completed.add(-1);
        setsDone = 0;
    }

    public void setSets(int sets) {
        this.sets = sets;
    }

    public void setReps(int reps) {
        this.reps = reps;
    }

    public void setStream(ArrayList<Float> accelerationStream){
        this.accelerationStream = accelerationStream;
    }

    public void setName(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
        return name + ": " + reps + "x" + sets;
    }

}
