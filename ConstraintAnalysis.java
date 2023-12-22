package falconsTools;

public class ConstraintAnalysis {
	static double StaticThrust = 49.45 ;
	static double CLmax = 1.8 ; 
	static double Wmax = 1 , Smin = 0 ;
	static double fCd , fCl , fSRef ;
	public static void main(String args[]) {
		double www  =fuselageSimulator ( ) ;
		findWS() ;		
		
	}

	public static void findWS( ) { 
		
		for ( double W = Wmax - 1 ; W <= 245 ; W += 1 ) { 
			for ( double S = Smin + 0.01 ; S <= 2 ; S += 0.01 ) { 
				double stallSpeed = stallSpeed( S , W ) ; 
				double takeOffSpeed = 1.1 * stallSpeed ;
				double T = Thrust ( takeOffSpeed) ; 
				double climbGradient = 0.1 ;
				double LDclimb = L_Dclimb(S) ;
				if ( equationSatisfied( T , W , LDclimb , climbGradient )) {
					System.out.println(S + "\t" + W ) ; 
					break ; }
					if ( W > Wmax  ) {
						System.out.println (Wmax) ;
						Smin = S ;
						Wmax = W ;
					
				}
			}
			
		}
		
		for ( double W = Wmax - 1 ; W <= 245 ; W += 1 ) { 
			for ( double S = 2 ; S >= Smin ; S -= 0.01 ) { 
				double stallSpeed = stallSpeed( S , W ) ; 
				double takeOffSpeed = 1.1 * stallSpeed ;
				double T = Thrust ( takeOffSpeed) ; 
				double climbGradient = 0.1 ;
				double LDclimb = L_Dclimb(S) ;
				if ( equationSatisfied( T , W , LDclimb , climbGradient )) {
					System.out.println(S + "\t" + W ) ; 
					break ; }
					if ( W > Wmax  ) {
						System.out.println (Wmax) ;
						Smin = S ;
						Wmax = W ;
					
				}
			}
			
		}
	
		double printParameters = optimiseSemiTaperedWing ( Smin , true ) ;	

		if ( printParameters == -1 ) { 
			System.out.println("S" + Smin);
			findWS() ; 
		}
		else {
		System.out.println("Maximum Take-Off Weight : " + Wmax/10 + " kg") ;
		System.out.println("Required Surface Area : " + Smin + " m^2") ;
		System.out.println("Stall Speed : " + stallSpeed(Smin, Wmax ) + " m/s") ;
		System.out.println ( "Thrust at take off : " + Thrust(1.1*stallSpeed(Smin,Wmax)) + " N") ;
		}
		
	}
	
	public static double Thrust ( double velocity ) {
		double T = StaticThrust  -0.0436348*velocity*velocity - 0.660454 * velocity ;
		return Math.max(T, 0);
	}
	
	public static double stallSpeed ( double S , double weight ) { 
		double stallSpeed  = 0 ;
		stallSpeed = Math.pow((weight /  ( 0.5 * 1.225 * S * CLmax ) ),0.5) ; 
		return stallSpeed ; 
	}
	public static double L_Dclimb ( double S ) {
		double L_Dclimb = optimiseSemiTaperedWing ( S , false ) ;
		return L_Dclimb ;
	}
	public static double fuselageSimulator ( ) { 
		double airfoilLiftSlope = 0.95, ZeroLiftAoA = -6.5, maxthickness = 20.5, maxthicknesspos = 29.2;

        // Flight Conditions
        double density = 1.225, Velocity = 13, Re, mu = 1.86 * 1 / 100000, Ma = Velocity / 343, pi = 3.14159265359;

        // Input
        double fuselageEfficiencyDecrease = 1; // set arbitrarily

        // Output
        double mc = 0, mb = 0, mlambda = 0, mphi = 0, maoi = 0, mLift = 0, mrange = 0, mendurance = 0, mscore = 0, mAR = 0,
                tipChord = 0, c = 0, b = 0, lambda = 0, phi = 0, aoiValue = 0, Lift = 0, range = 0, endurance = 0,
                Sref = 0, sweepAngle = 0, AR = 0, Ne = 0, TR = 0, flambda = 0, e = 0, Cf = 0, FF = 0, Swet = 0, Cd0 = 0,
                liftSlope = 0, yIntercept = 0, ceff = 0, S1 = 0, S2 = 0, Cdi = 0, Cd = 0, CLclimb = 0, score = 0,
                msweepAngle = 0, mSRef = 0,myIntercept = 0, mliftSlope = 0;
        aoiValue = 8 ;
        ceff = 0.86 ;
        b = 3 ;
        Sref = ceff * 0.4 ;
        AR = ceff / b ;
		Re = ceff * density * Velocity / mu;        
        TR = 1 ; 
        flambda = 0.005 * (1 + 1.5 * Math.pow((TR - 0.6), 2));
        sweepAngle =0;
        e = 0.8 ;

        Cf = 0.074 / Math.pow(Re, 0.2) ;
        FF = ((1 + (maxthicknesspos / 100 > 0.3 ? 1.2 : 2) * maxthickness / 100 +
                100 * Math.pow(maxthickness / 100, 4)) * 1.05) ;
        Swet = 2 * (1 + 0.25 * maxthickness / 100) ;
        Cd0 = Swet * FF * Cf;
        liftSlope = airfoilLiftSlope / (1 + 57.3 * airfoilLiftSlope / (pi * e * AR));
        yIntercept = -ZeroLiftAoA * liftSlope;
        CLclimb = ( liftSlope * aoiValue + yIntercept ) * 0.8;
        Cdi = CLclimb * CLclimb / (pi * e * AR);
        Cd = Cd0 + Cdi;
        Lift = 0.5 * density * Velocity * Velocity * Sref * CLclimb;
        range = CLclimb / Cd;
        endurance = Math.pow(CLclimb, 1.5) / Cd;
        fCd = Cd ;
        fCl = CLclimb ;
        fSRef = Sref ;
        return Cd*0.5 * 1.225 * Velocity * Velocity * Sref * 1;
	}
	
    public static double optimiseSemiTaperedWing(double S , boolean printParameters ) {
        // 2D Airfoil Data
        double airfoilLiftSlope = 0.10505, ZeroLiftAoA = -14.2, maxthickness = 14.53, maxthicknesspos = 25.14;

        // Flight Conditions
        double density = 1.225, Velocity = 13, Re, mu = 1.86 * 1 / 100000, Ma = Velocity / 343, pi = 3.14159265359;

        // Input
        double fuselageEfficiencyDecrease = 1; // set arbitrarily

        // Output
        double mc = 0, mb = 0, mlambda = 0, mphi = 0, maoi = 0, mLift = 0, mrange = 0, mendurance = 0, mscore = 0, mAR = 0,
                tipChord = 0, c = 0, b = 0, lambda = 0, phi = 0, aoiValue = 0, Lift = 0, range = 0, endurance = 0,
                Sref = 0, sweepAngle = 0, AR = 0, Ne = 0, TR = 0, flambda = 0, e = 0, Cf = 0, FF = 0, Swet = 0, Cd0 = 0,
                liftSlope = 0, yIntercept = 0, ceff = 0, S1 = 0, S2 = 0, Cdi = 0, Cd = 0, CLclimb = 0, score = 0,
                msweepAngle = 0, mSRef = 0,myIntercept = 0, mliftSlope = 0 ,mdrag = 0 , me = 0;

        for (aoiValue = 10; aoiValue <= 10; aoiValue += 1) {
            for (c = 0.4; c <= 0.8; c += 0.005) {
                for (b = 2.4; b <= 2.6; b += 0.01) {
                    for (lambda = 1; lambda <= 1; lambda += 0.05) {
                        for (phi = 1; phi <= 1; phi += 0.05) {
                            tipChord = lambda * c ;
                            if (tipChord < 0.36)
                                break;
                            S1 = b * c * phi;
                            S2 = b * c * (1 + lambda) / 2 * (1 - phi);
                            Sref = S1 + S2;
                            ceff = Sref / b;
                            Re = ceff * density * Velocity / mu;
                            AR = b * b / Sref;
                            TR = (S1 + S2 * lambda) / Sref;
                            flambda = 0.005 * (1 + 1.5 * Math.pow((TR - 0.6), 2));
                            sweepAngle = Math.atan((-0.25 * c + 0.25 * lambda * c) / (b / 2));
                            e = 1 / ((1 + 0.12 * Ma * Ma) *
                                    ((1 + ((0.142 + flambda * AR * Math.pow(10 * maxthickness / 100, 0.33)) /
                                            Math.pow(Math.cos(sweepAngle), 2)) +
                                            (0.3 * Ne + 0.1) / Math.pow(4 + AR, 0.8))));

                            Cf = 0.074 / Math.pow(Re, 0.2);
                            FF = ((1 + (maxthicknesspos / 100 > 0.3 ? 1.2 : 2) * maxthickness / 100 +
                                    100 * Math.pow(maxthickness / 100, 4)) * 1.05);
                            Swet = 2 * (1 + 0.25 * maxthickness / 100);
                            Cd0 = Swet * FF * Cf;
                            liftSlope = airfoilLiftSlope / (1 + 57.3 * airfoilLiftSlope / (pi * e * AR));
                            yIntercept = -ZeroLiftAoA * liftSlope;
                            CLclimb = ( liftSlope * aoiValue + yIntercept ) * 0.8;
                            CLclimb = (CLclimb * Sref + fCl * fSRef ) / ( Sref + fSRef ) ;
                            Cdi = CLclimb * CLclimb / (pi * e * AR);
                            Cd = (Cd0 + Cdi);
                            Cd = (Cd * Sref + fCd * fSRef ) / ( Sref + fSRef ) ;
                            Lift = 0.5 * density * Velocity * Velocity * Sref * CLclimb;
                            range = CLclimb / Cd;
                            endurance = Math.pow(CLclimb, 1.5) / Cd;

                            score = range * 0.4 + endurance * 0.6 - fuselageEfficiencyDecrease * c;
                            if (Sref >= S ) {

                                if (range > (mrange)) {
                                	me =e ;
                                    mb = b;
                                    mc = c;
                                    mphi = phi;
                                    mlambda = lambda;
                                    mendurance = endurance;
                                    mrange = range;
                                    maoi = aoiValue;
                                    mLift = Lift;
                                    mscore = score;
                                    mAR = AR;
                                    msweepAngle = sweepAngle;
                                    mSRef = Sref;
                                    mliftSlope = liftSlope ; 
                                    myIntercept = yIntercept ; 
                                    mdrag = Cd0 ;
                                }
                            }
                        }
                    }
                }
            }
        }
        if ( printParameters ) { 
        	double Creality = 0.8 ;
        	double TaildownforceFactor = 0.95 ;
        	
        	double CLcruise = ( mliftSlope * 5 + myIntercept ) * Creality * TaildownforceFactor ;
    	    double Vcruise =  Math.pow(( Wmax /  ( 0.5 * 1.225 * mSRef * CLcruise ) ),0.5)  ;
    	    double CDcruise = 0.5 * 1.225 * Vcruise * Vcruise * mSRef * ( mdrag + CLcruise * CLcruise / (pi * me * mAR)) ;
    	    //if ( Vcruise > ( 1.3 * stallSpeed ( Smin , Wmax)) )
        	   // return -1 ;
        	System.out.println("Wingspan = " + String.format("%.3f", mb));
    	    System.out.println("Rootchord = " + String.format("%.3f", mc));
    	    System.out.println("Taper Start Point = " + String.format("%.3f", mphi));
    	    System.out.println("Taper Ratio = " + String.format("%.3f", mlambda));
    	    System.out.println("Tip chord = " + String.format("%.3f", mc*mlambda));
    	    System.out.println("mendurance = " + String.format("%.3f", mendurance));
    	    System.out.println("mrange = " + String.format("%.3f", mrange));
    	    System.out.println("maoi = " + String.format("%.3f", maoi));
    	    System.out.println("mLift = " + String.format("%.3f", mLift));
    	    System.out.println("AR = " + String.format("%.3f", mAR));
    	    System.out.println("wing Sweep = " + String.format("%.3f", msweepAngle*57.3)+ " deg");
    	    System.out.println ( "SRef : " + mSRef) ;
    	    System.out.println("Drag afforded :" + CDcruise) ;
    	    
    	    
    	    System.out.println ( "cruise velocity : " + Vcruise + " m/s") ;
    	    System.out.println ( "Sref: " + mSRef ) ;
    	    System.out.println ( "Wmax : " + Wmax ) ;
    	    System.out.println ( "CLcruise : " + CLcruise ) ;
    	    System.out.println ( fuselageSimulator () ) ;
    	    
    	    return 1 ;
        }
        return mrange ;
    }
    
    public static boolean equationSatisfied ( double T , double W , double LDclimb , double climbGradient ) {
    	climbGradient = 0.15 ;
    	return ( ( T/W ) >= (1 / (LDclimb) + climbGradient ));
    }
}
