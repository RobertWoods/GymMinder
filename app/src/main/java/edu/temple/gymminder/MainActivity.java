package edu.temple.gymminder;

import android.content.Intent;
import android.os.Bundle;
import android.support.annotation.NonNull;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v4.app.FragmentTransaction;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;

import com.google.firebase.auth.FirebaseAuth;
import com.google.firebase.auth.FirebaseUser;

import edu.temple.gymminder.geofence.GeofenceFragment;

public class MainActivity extends AppCompatActivity implements SigninFragment.SigninListener,
        MainFragment.DetailListener, WorkoutCreatorFragment.Listener, AdHocCreatorFragment.Listener{

    public static final String AD_HOC = "Laughing to the bank like ahhHA";
    public static final String START_FRAGMENT_EXTRA = "It was always me vs the world." +
            "Until I found it was me vs me.";

    private FirebaseAuth auth;
    private Fragment activeFragment;
    
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        setupAuth();
        if (BuildConfig.FLAVOR.equals("espresso")) {
            auth.signOut();
        }
        FragmentManager fragmentManager = getSupportFragmentManager();
        if (auth.getCurrentUser() == null) {
            fragmentManager.beginTransaction()
                    .replace(R.id.mainFrame, new SigninFragment())
                    .commit();
        } else if(getIntent().getExtras()!=null) {
            if(getIntent().getExtras().get(START_FRAGMENT_EXTRA)!=null){
                handleStartFragmentExtra(getIntent().getExtras());
            }
        } else {
            goToMain();
        }
    }

    @Override
    protected void onNewIntent(Intent intent) {
        Bundle extras = intent.getExtras();
        if(extras!=null) {
            if(extras.get(START_FRAGMENT_EXTRA)!=null) {
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
                startFragment(new MainFragment());
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
                break;
//            case R.id.geofenceOption:
//                startFragment(new GeofenceFragment());
//                break;
        }
        return true;
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
                    getSupportFragmentManager().beginTransaction().replace(R.id.mainFrame,
                            new SigninFragment()).commit();
                }
            }
        });
    }

    public void goToMain() {
        startFragment(new MainFragment());
    }

    public void goToDetail(Workout workout, String name) {
        DetailFragment detailFragment = DetailFragment.newInstance(workout, name);
        startFragment(detailFragment);
    }

    public void goToWorkoutCreator() {
        WorkoutCreatorFragment workoutCreatorFragment = new WorkoutCreatorFragment();
        getSupportFragmentManager().beginTransaction()
                .replace(R.id.mainFrame, workoutCreatorFragment)
                .addToBackStack(null)
                .commit();
    }

    @Override
    public void finishFragment(Fragment f) {
        getSupportFragmentManager().popBackStack();
    }


    public void startFragment(Fragment fragment){
        FragmentTransaction transaction = getSupportFragmentManager()
                .beginTransaction()
                .replace(R.id.mainFrame, fragment);
        if(activeFragment instanceof MainFragment) transaction = transaction.addToBackStack(null);
        transaction.commit();
        activeFragment = fragment;
    }


}
