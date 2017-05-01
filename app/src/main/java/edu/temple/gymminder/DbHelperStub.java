package edu.temple.gymminder;

import com.google.firebase.auth.FirebaseUser;

import java.util.ArrayList;
import java.util.Date;

import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;

/**
 * Created by rober_000 on 4/11/2017.
 */

public class DbHelperStub extends DbHelper {


    public DbHelperStub(Listener listener) {
        super(listener);
    }

    @Override
    public void retrieveAllWorkouts(FirebaseUser user) {
        ArrayList<Workout> workouts = new ArrayList<>();
        ArrayList<Exercise> exercises = new ArrayList<>();
        exercises.add(new Exercise("Bench", 3, 3));
        workouts.add(new Workout(exercises));
        ArrayList<String> names = new ArrayList<>();
        names.add("hey p3p");
        if(listener!=null) {
            listener.respondToWorkouts(workouts, names);
        }
    }

    @Override
    public void getCatalog() {
        ArrayList<Exercise> exercises = new ArrayList<>();
        exercises.add(new Exercise("Curls", -1, -1));
        exercises.add(new Exercise("Bench", -1, -1));
        exercises.add(new Exercise("Squat", -1, -1));
        if (listener!=null) {
            listener.respondToCatalog(exercises);
        }
    }

    @Override
    public void addWorkout(Workout workout, String workoutName, FirebaseUser user, Date date){
        if(listener!=null) {
            listener.onWorkoutAdded();
        }
    }

    @Override
    public void addNewWorkout(Workout workout, String workoutName, FirebaseUser user){
        if(listener!=null) {
            listener.onWorkoutAdded();
        }
    }



}
