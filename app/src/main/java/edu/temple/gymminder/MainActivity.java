package edu.temple.gymminder;

import android.content.Intent;
import android.os.Bundle;
import android.support.annotation.NonNull;
import android.support.design.widget.BottomNavigationView;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v4.app.FragmentTransaction;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.view.View;

import com.google.firebase.auth.FirebaseAuth;
import com.google.firebase.auth.FirebaseUser;

import java.util.ArrayList;

import java.util.Date;

import edu.temple.gymminder.geofence.GeofenceFragment;
import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;

public class MainActivity extends AppCompatActivity
        implements SigninFragment.SigninListener,
        WorkoutsFragment.DetailListener,
        WorkoutCreatorFragment.Listener,
        AdHocCreatorFragment.Listener,
        HistoryFragment.OnFragmentInteractionListener {

    public static final String AD_HOC = "Laughing to the bank like ahhHA";
    public static final String START_FRAGMENT_EXTRA = "It was always me vs the world." +
            "Until I found it was me vs me.";
    private static final String MAIN_FRAGMENT_TAG = "System system blower";
    private static final String OTHER_FRAGMENT_TAG = "Grab your fucking chain and drag it through the bike lane";

    private FirebaseAuth auth;
    private Fragment activeFragment;
    private BottomNavigationView tabBar;

    
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        setupAuth();

        tabBar = (BottomNavigationView) findViewById(R.id.navigation);
        tabBar.setOnNavigationItemSelectedListener(new BottomNavigationView.OnNavigationItemSelectedListener() {
            @Override
            public boolean onNavigationItemSelected(@NonNull MenuItem item) {
                switch (item.getItemId()) {
                    case R.id.navigation_workouts:
                        goToWorkouts();
                        return true;
                    case R.id.navigation_history:
                        goToHistory();
                        return true;
                    case R.id.navigation_geofence:
                        goToGeofence();
                        return true;
                    default:
                        return false;
                }
            }
        });

        if (BuildConfig.FLAVOR.equals("espresso")) {
            auth.signOut();
            DataUtils.loadRepetitionFile("Bench", this).delete();
            DataUtils.loadRepetitionFile("Deadlift", this).delete();
            DataUtils.loadRepetitionFile("Squat", this).delete();
            DataUtils.loadRepetitionFile("Curl", this).delete();
        }
        if (auth.getCurrentUser() == null) {
            startFragment(new SigninFragment());
        } else if (getIntent().getExtras()!=null) {
            if (getIntent().getExtras().get(START_FRAGMENT_EXTRA) != null){
                handleStartFragmentExtra(getIntent().getExtras());
            } else {
                goToWorkouts();
            }
        } else {
            goToWorkouts();
        }
    }

    @Override
    protected void onNewIntent(Intent intent) {
        Bundle extras = intent.getExtras();
        if(extras!=null) {
            if (extras.get(START_FRAGMENT_EXTRA)!=null) {
                handleStartFragmentExtra(extras);
            }
        }
    }

    private void handleStartFragmentExtra(Bundle extras){
        String fragment = extras.getString(START_FRAGMENT_EXTRA, "");
        switch (fragment) {
            case AD_HOC:
                startFragment(new AdHocCreatorFragment());
                break;
            default:
                startFragment(new WorkoutsFragment());
        }
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        MenuInflater inflater = getMenuInflater();
        inflater.inflate(R.menu.options_menu, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        switch (item.getItemId()){
            case R.id.signOutOption:
                auth.signOut();
                tabBar.setVisibility(View.GONE);
                break;
        }
        return true;
    }

    @Override
    public void onBackPressed(){
        super.onBackPressed();
        FragmentManager fragmentManager = getSupportFragmentManager();
        activeFragment = fragmentManager.findFragmentByTag(MAIN_FRAGMENT_TAG);
        tabBar.getMenu().findItem(R.id.navigation_workouts).setChecked(true);
        /*
            If our fragment depth (number of fragments navigating away from main fragment)
            becomes higher than 2, we want to remove any fragments not handled by the back stack.
         */
        Fragment oldFragment = fragmentManager.findFragmentByTag(OTHER_FRAGMENT_TAG);
        if(oldFragment!=null) {
            fragmentManager
                    .beginTransaction()
                    .remove(oldFragment)
                    .commit();
        }
    }

    public void setupAuth() {
        auth = FirebaseAuth.getInstance();
        auth.addAuthStateListener(new FirebaseAuth.AuthStateListener() {
            @Override
            public void onAuthStateChanged(@NonNull FirebaseAuth firebaseAuth) {
                FirebaseUser user = firebaseAuth.getCurrentUser();
                if (user != null) {
                    Log.d("Auth", "Signed in");
                } else {
                    Log.d("Auth", "Not Signed In");
                    startFragment(new SigninFragment());
                }
            }
        });
    }

    public void goToWorkouts() {
        startFragment(new WorkoutsFragment());
        tabBar.setVisibility(View.VISIBLE);
        tabBar.getMenu().findItem(R.id.navigation_workouts).setChecked(true);
    }

    public void goToWorkoutsDetail(Workout workout, String name) {
        DetailFragment detailFragment = DetailFragment.newInstance(workout, name);
        startFragment(detailFragment);
    }

    public void goToWorkoutHistoryDay(Date day){
        //DayFragment dayFragment = DayFragment.newInstance(day);

    }

    public void goToWorkoutCreator() {
        WorkoutCreatorFragment workoutCreatorFragment = new WorkoutCreatorFragment();
        startFragment(workoutCreatorFragment);
    }

    @Override
    public void goToAdHocCreator(ArrayList<Exercise> exercises, int selected) {
        AdHocCreatorFragment fragment = AdHocCreatorFragment.newInstance(exercises, selected);
        startFragment(fragment);
    }

    public void goToHistory() {
        startFragment(new HistoryFragment());
        tabBar.setVisibility(View.VISIBLE);
        tabBar.getMenu().findItem(R.id.navigation_history).setChecked(true);
    }

    public void goToGeofence() {
        startFragment(new GeofenceFragment());
        tabBar.setVisibility(View.VISIBLE);
        tabBar.getMenu().findItem(R.id.navigation_geofence).setChecked(true);
    }


    @Override
    public void finishFragment(Fragment f) {
        getSupportFragmentManager().popBackStack();
    }


    public void startFragment(Fragment fragment){
        String tag = fragment instanceof WorkoutsFragment ? MAIN_FRAGMENT_TAG : OTHER_FRAGMENT_TAG;
        FragmentTransaction transaction = getSupportFragmentManager()
                .beginTransaction()
                .replace(R.id.mainFrame, fragment, tag);
        if(activeFragment instanceof WorkoutsFragment) transaction = transaction.addToBackStack(null);
        transaction.commit();
        activeFragment = fragment;
    }

}
