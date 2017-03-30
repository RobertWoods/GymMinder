package edu.temple.gymminder;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Created by rober_000 on 1/31/2017.
 */

public class Workout implements Serializable {
    ArrayList<Exercise> exercises;

    public Workout() {

    }

    public Workout(ArrayList<Exercise> exercises) {
        this.exercises = exercises;
    }

    public Workout(ArrayList<Exercise> exercises, boolean errorCheck) {
        for (Exercise e : exercises) {
            if (e.workout.length() > 50 || e.reps <= 0 || e.sets <= 0) exercises.remove(e);
        }
        this.exercises = exercises;
    }

    @Override
    public String toString() {
        String res = "";
        for (Exercise e : exercises) {
            res += e.workout + ": " + e.reps + "x" + e.sets + "\n";
        }
        res = res.substring(0, res.length() - 1);
        return res;
    }


}
