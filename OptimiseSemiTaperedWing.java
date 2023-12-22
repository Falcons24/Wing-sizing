package falconsTools;
public class OptimiseSemiTaperedWing {
	
	//2D Airfoil Data
	static double airfoilLiftSlope = 0.10505 , ZeroLiftAoA = -14.2 , maxthickness = 14.53 , maxthicknesspos = 25.14 ;
	
	//Flight COnditions
	static double density = 1.225 , Velocity = 13 , Re , mu = 1.86 * 1/100000 , Ma = Velocity / 343 , pi = 3.14159265359;
	
	//Input
	static double targetLift = 210 ;
			; 
	static double fuselageEfficiencyDecrease = 1 ; //set arbitrarily
	
	//Output 
	static double mc , mb , mlambda , mphi , maoi , mLift , mrange , mendurance ,mscore , mAR , tipChord ; 
	
	//variables
	static double c , b , lambda , phi , aoi ,Lift , range , endurance ; 
	
	//Intermediate values
	static double Sref , sweepAngle = 0, AR, Ne =0, TR, flambda, e , Cf, FF, Swet, Cd0, liftSlope , yIntercept , ceff, S1 , S2,
			Cdi ,Cd , CLcruise , score , msweepAngle , mSRef; 
	
	public static void main ( String args[] ) { 
		//LiftCalculator ( 0.6,2.6,5,0.24,0.46 );
		for ( aoi = 5 ; aoi <= 5 ; aoi += 1) {
		for ( c = 0.6; c <= 0.65 ; c += 0.005 ) {
			for (  b = 2.6; b <= 2.6 ; b += 0.01 ) {
				for (  lambda = 0 ; lambda <= 0.7; lambda+= 0.01) { 
					for ( phi = 0.3;  phi <= 0.5 ; phi+= 0.01 ) { 
						tipChord = lambda*c ;
						if ( tipChord < 0.36 ) 
							break ;
						S1 = b * c * phi ;
						S2 = b * c * ( 1 + lambda ) /2 * ( 1 - phi ) ;
						Sref = S1 + S2 ; 
						ceff = Sref / b ;
						Re = ceff * density * Velocity / mu ; 
						AR = b*b / Sref ; 
						TR = ( S1 + S2*lambda ) / Sref ; 
						flambda =  0.005* ( 1 + 1.5 * Math.pow( (TR - 0.6) , 2 ) ) ;
						sweepAngle = Math.atan((-0.25 * c + 0.25 * lambda * c) / (b / 2));
						e = 1/((1 + 0.12 * Ma*Ma) *
				                ((1 + ((0.142 + flambda * AR * Math.pow(10 * maxthickness / 100, 0.33)) / Math.pow(Math.cos(sweepAngle), 2)) +
				                (0.3 * Ne + 0.1) / Math.pow(4 + AR, 0.8))));
								
						Cf = 0.074 / Math.pow( Re,0.2 ) ;
						FF = ((1 + (maxthicknesspos / 100 > 0.3 ? 1.2 : 2) * maxthickness / 100 + 100 * Math.pow(maxthickness / 100, 4)) * 1.05);
						Swet = 2*(1+0.25* maxthickness / 100 ) ;
						Cd0 = Swet * FF * Cf ;
						liftSlope =airfoilLiftSlope/(1+57.3*airfoilLiftSlope/(pi*e*AR) );
						yIntercept= -ZeroLiftAoA * liftSlope ; 
						CLcruise = liftSlope * aoi+yIntercept ; 
						Cdi = CLcruise * CLcruise / (pi*e*AR) ;
						Cd = Cd0 + Cdi ;
						Lift = 0.5 * density * Velocity * Velocity * Sref * CLcruise ; 
						range = CLcruise / Cd ; 
						endurance = Math.pow ( CLcruise , 1.5 )/ Cd ; 	
						
						
						score = range*0.4+endurance*06 - fuselageEfficiencyDecrease*c  ;
						if ( Sref >= 1.51 ) { 
							
							if ( score > ( mscore ) ) { 
								mb = b ; 
								mc = c ; 
								mphi = phi ; 
								mlambda = lambda ; 
								mendurance = endurance ; 
								mrange = range ;
								maoi = aoi ;
								mLift = Lift ;								
								mscore = score ;
								mAR = AR ;
								msweepAngle = sweepAngle ;
								mSRef = Sref ;
							}
						}
					}
				}
			}
		}
	}
		output() ;
	}
	public static void LiftCalculator ( double c , double b , double aoi , double phi , double lambda ) {
		S1 = b * c * phi ;
		S2 = b * c * ( 1 + lambda ) /2 * ( 1 - phi ) ;
		Sref = S1 + S2 ; 
		ceff = Sref / b ;
		Re = ceff * density * Velocity / mu ; 
		AR = b * b / Sref ; 
		TR = ( S1 + S2 * lambda ) / Sref ; 
		flambda =  0.005* ( 1 + 1.5 * Math.pow ( (TR - 0.6) , 2 ) ) ;
		sweepAngle = Math.atan((-0.25 * c + 0.25 * lambda * c) / (b / 2));
		e = 1/((1 + 0.12 * Ma* Ma) *
                ((1 + (0.142 + flambda * AR * Math.pow(10 * maxthickness / 100, 0.33)) / Math.pow(Math.cos(sweepAngle), 2) +
                (0.3 * Ne + 0.1) / Math.pow(4 + AR, 0.8))));
		
		Cf = 0.074 / Math.pow( Re,0.2 ) ;
		FF = ((1 + (maxthicknesspos / 100 > 0.3 ? 1.2 : 2) * maxthickness / 100 + 100 * Math.pow(maxthickness / 100, 4)) * 1.05);
		Swet = 2*(1+0.25*maxthickness/100) ;
		Cd0 = Swet * FF * Cf ;
		liftSlope =airfoilLiftSlope/(1+57.3*airfoilLiftSlope/(pi*e*AR) );
		yIntercept=-ZeroLiftAoA*liftSlope ; 
		CLcruise = liftSlope*aoi+yIntercept ; 
		Cdi =CLcruise*CLcruise/(pi*e*AR) ;
		Cd = Cd0 + Cdi ;
		Lift = 0.5 * density * Velocity * Velocity * Sref * CLcruise ; 
		range = CLcruise / Cd ; 
		endurance = Math.pow ( CLcruise , 1.5 )/ Cd ; 
		
		System.out.println("S1 = " + S1);
		System.out.println("S2 = " + S2);
		System.out.println("Sref = " + Sref);
		System.out.println("ceff = " + ceff);
		System.out.println("Re = " + Re);
		System.out.println("AR = " + AR);
		System.out.println("TR = " + TR);
		System.out.println("flambda = " + flambda);
		System.out.println("sweepAngle = " + sweepAngle);
		System.out.println("e = " + e);
		System.out.println("Cf = " + Cf);
		System.out.println("FF = " + FF);
		System.out.println("Swet = " + Swet);
		System.out.println("Cd0 = " + Cd0);
		System.out.println("liftSlope = " + liftSlope);
		System.out.println("yIntercept = " + yIntercept);
		System.out.println("CLcruise = " + CLcruise);
		System.out.println("Cdi = " + Cdi);
		System.out.println("Cd = " + Cd);
		System.out.println("Lift = " + Lift);
		System.out.println("range = " + range);
		System.out.println("endurance = " + endurance);
	}	
	public static void output() {
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
	    System.out.println("MGC : "+ mSRef / mb) ;
    }
}
