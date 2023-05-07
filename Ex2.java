

import java.util.Arrays;



/**
 * Introduction to Computer Science 2023, Ariel University,
 * Ex2: arrays, static functions and JUnit
 *
 * This class represents a set of functions on a polynom - represented as array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynom: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code here")
 *
 * @author boaz.benmoshe
 */
public class Ex2 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynom is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynom at x.
	 * @param poly
	 * @param x
	 * @return f(x) - the polynom value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans +=c*poly[i];
		}
		return ans;
	}

	/** Given a polynom (p), a range [x1,x2] and an epsilon eps. 
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x1) <= 0. 
	 * This function should be implemented recursively.
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double f2 = f(p,x2);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (f1*f2<=0 && Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 * ------------
	 * the code is based on the linear equation: y= mx+b---> in our condtion y=bx+c.
	 * - in the case we have length equal to 3 (we have 3 points)
	 * the code is based in quadratic equation: y=ax^2+bx+c
	 * in each case we define a new array of ans. 
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length; //the length of xx array
		int ly = yy.length; //the length of yy array

		//lx should be equal to ly and or 2 or 3
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
			// check the case of two points
			//(x1,y1),(x2,y2)
			if (lx==2 ) {
				double x1=xx[1]	,x2=xx[0];
				double y1=yy[1],y2=yy[0];
				//linear equation
				double a= 0;
				double b=(y2 - y1) / (x2 - x1);
				double c = y1 - b * x1;
				ans= new double [2];
				ans[0]=b;
				ans[1]=c;

			}
			//check the case of three points
			//(x1,y1),(x2,y2),(x3,y3)
			else if (lx==3) {
				double x1=xx[2]	,x2=xx[1],x3=xx[0];
				double y1=yy[2],y2=yy[1],y3=yy[0];

				//quadratic equation
				double	divider = (x1-x2) * (x1-x3) * (x2-x3);
				double	a     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / divider;
				double	b     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / divider;
				double	c    = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / divider;

				//ans = new double[]{c,b,a};
				ans = new double[3];
				ans[0] = c;
				ans[1] = b;
				ans[2] = a;
			}
		}

		return ans;
	}

	/** Two polynoms are equal if and only if the have the same values f(x) for 1+n values of x, 
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynom
	 * @param p2 second polynom
	 * @return true iff p1 represents the same polynom as p2.
	 * --------------------------------
	 * explantion:
	 *  in this function first we Compute the maximum degree n of the two polynoms. 
	 * then Evaluate the two polynomsat n+1 values of  variable x then we set x to be i
	 * using "f function" to check if value of the the two polynoms at x are equal.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans=true;	
		int n= Math.max(p1.length, p2.length);
		for (int i=0;i<=n;i++) {//iterate through the values to n+1 
			int x=i;//initialize x to be i
			if (Math.abs(f(p1,x)-f(p2,x))>=EPS) {//using f function and check the case
				return false;
			}

		}
		return ans;

	}


	/** 
	 * Computes a String representing the polynom.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynom represented as an array of doubles
	 * @return String representing the polynom: 
	 */

	public static String poly(double[] poly) {
		String ans = "";//initialize an empty string
		if(poly.length==0) {//check if the 'poly' is empty
			ans="0";
		}
		else if(poly.length==1) {//check if the length of 'poly'is 1
			ans+=poly[0];//set the single  coefficient 
		}
		else {// if the length of 'poly' greater than 1 then iterate through its coefficients.
			for(int i= poly.length-1; i>=0; i--) {
				if(poly[i]!=0) {

					if(ans.isEmpty()) {// Check if the string ans is empty to deal with the signal ,if it is we dont need to set - or +.
						if(poly[i]<0) {
							ans += "-";
						}
					}
					else {
						if(poly[i]>0) {//check if the current  coefficient is positive to set plus'+'signal to ans.
							ans += " + ";
						}
						else if (poly[i]<0){//check if the current  coefficient is negative to set minus'-'signal to ans.
							ans += " - ";
						}
					}
					if(Math.abs(poly[i])!=1||i==0) {//Check if the absolute value of the current coefficient is not 1 or in the index 0
						ans += Math.abs(poly[i]);
					}
					if(i>0) {//Check if the current term has an x variable
						ans += "x";
					}
					if(i>1) {//Check if the current term has an exponent greater than 1 to set the power
						ans += "^" + i;
					}
				}
			}
		}
		return ans;
	}

	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {

		double y1_1 = f(p1,x1) , y1_2 = f(p1,x2); //the y-values of the first polynomial function at the endpoints x1 and x2
		double y2_1 = f(p2,x1) , y2_2 = f(p2,x2);//the y-values of the second polynomial function at the endpoints x1 and x2

		//find the midpoint of the interval x1 and x2
		//and find the y-value of the first polynomial function at the midpoint mid
		double mid = (x1+x2)/2 , y1_12 = f(p1,mid);

		//the y-value of the second polynomial function at the midpoint mid
		double y2_12 = f(p2,mid);

		// here we find	the difference between the y-values of 
		//the first and second polynomial functions at x1 and x2

		double del3 = (y1_12 - y2_12);//the difference between the y-values of the first and second polynomial functions at the midpoint mid
		double del1 = (y1_1 - y2_1), del2 = (y1_2 - y2_2);
		while ((del1*del2<=0 )&& (Math.abs(del3)>=EPS)) {
			if(del3*del1<=0) {
				x2=mid;// define x2 to be mid
				y1_2= f(p1,x2);
				y2_2= f(p2,x2);
			}
			else {
				x1=mid;// define x1 to be mid
				y1_1 = f(p1,x1);y2_1 = f(p2,x1);
			}
			mid = (x1+x2)/2;
			del1 = (y1_1 - y2_1);
			del2 = (y1_2 - y2_2);
			y1_12 = f(p1,mid);
			y2_12 = f(p2,mid);
			del3 = (y1_12 - y2_12);
		}
		return mid;
	}	

	/**
	 * Given a polynom (p), a range [x1,x2] and an integer with the number (n) of sample points. 
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 * -------------------------------
	 * in this code we deal with find the lingth by considering each segment is side of triangle 
	 * and each time we find the thre sides of triangle 
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
		double seg =Math.abs(x1-x2)/numberOfSegments;// calculate the segment length
		double [] a =new double[numberOfSegments+1];// initialize array for widths values like x
		double []b=new double[numberOfSegments+1];// initialize array for heights values like y


		// calculate a and b values for each sample point
		for(int i=0;i<=numberOfSegments;i++) {
			a[i] =x1+i*seg;
			b[i]=f(p,a[i]);
		}
		// calculate the length of each segment and add to the answer
		for(int i=0;i<numberOfSegments;i++) {
			//pythagoras equation a^2+b^2=c^2
			double aa=Math.pow(a[i+1]-a[i],2);
			double bb=Math.pow(b[i+1]-b[i],2);

			double cc=Math.sqrt(aa+bb);
			ans+=cc;//
		}
		return ans;
	}


	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom). 
	 * This function computes an approximation of the area between the polynoms within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynoms within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
		if(numberOfTrapezoid>=1) {//numberOfTrapezoid should be greater than 1
			double len=(x2-x1)/numberOfTrapezoid;//find the width of each Trapezoid
			for(double i=x1;i<x2;i+=len) {
				double len2=Math.abs(f(p1,i)-f(p2,i));//calulate the length in each value in  [x1,x2]
				ans+=len*len2;

			}
		}
		return ans;
	}


	/**
	 * This function computes the array representation of a polynom from a String
	 * representation. Note:given a polynom represented as a double array,  
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynom.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		if ((p != null) && (p != "")) {// Check if the input string is not null or empty 

			p = p.replaceAll(" ", "");//remove the spaces
			p = p.replaceAll("-", "+-");

			if (p.startsWith("+-")) {
				p = p.substring(1);
			}
			String[] p2 = p.split("\\+");//new string without the "+"
			int length;
			// find the degree of the polynomial based on the highest power of x in the terms
			if (p.contains("^")) {
				length = Integer.parseInt(p2[0].substring(p2[0].indexOf("^") + 1)) + 1;
			}
			else if (p.contains("x")) {//like 3x+1
				length =2;
			} else {//any number without variables like f(x)=3
				length = 1;
			}
			double[] ans = new double[length];
			int power;
			double coeff;
			// Iterate through the polynomial terms, extracting
			//the coefficient and power for each
			for (int i = 0; i < p2.length; i++) {
				if (p2[i].contains("x")) {// Term contains an x variable
					if (p2[i].contains("^")) {// Term contains a power
						power = Integer.parseInt(p2[i].substring(p2[i].indexOf("^") + 1));
						coeff = PolyCoefficients(p2[i]);//Extract the coefficient
					}
					else {
						power = 1;
						coeff =PolyCoefficients(p2[i]);//Extract the coefficient
					}
				} else {
					coeff = Double.parseDouble(p2[i]);//Extract the coefficient
					power = 0;
				}
				// Store the coefficient at the appropriate index in the coefficient array		
				ans[power] = coeff;
			}
			return ans;
		}
		return null;
	}
	/*this function will help us to do the getpolynomfromstring function
	 * in this function we get the coeffecitns of polynom represnts like String
	 * 	
	 */

	public static double PolyCoefficients(String poly) {
		if (poly.startsWith("x")) {// If the polynomial starts with "x", its coefficient is 1.
			return 1;	
		}
		else if (poly.startsWith("-x")) {// If the polynomial starts with "x", its coefficient is -1.
			return -1;	
		}
		else {// Otherwise, the coefficient is extracted from the substring up to the first occurrence of "x".

			return Double.parseDouble(poly.substring(0,poly.indexOf("x")));
		}
	}

	/**
	 * This function computes the polynom which is the sum of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return

	 */
	public static double[] add(double[] p1, double[] p2) {
		int k=p1.length,l=p2.length;
		int maxLength=Math.max(k,l);// Find the maximum length between the two arrays
		double [] ans = new double [maxLength];//to store the sum of the two polynomials

		// Loop over the arrays until the maximum length
		for(int i=0;i<maxLength;i++) {

			// Get the coefficient of p1 at index i
			double coefficients1=0.0; 
			if (i<p1.length) {
				coefficients1=p1[i];	
			}
			// Get the coefficient of p2 at index i
			double coefficients2=0.0;
			if(i<p2.length) {
				coefficients2=p2[i];	
			}
			// Compute the sum of the two coefficients and store it in ans at index i
			ans[i]=coefficients1+coefficients2;

		}
		ans = erasezero(ans);

		return ans;
	}
	/*
	 * this function help us in add function to remove
	 *  the zeros in the last of the array of polynom .
	 */
	private static double [] erasezero(double []pp) {
		while(pp[pp.length - 1] == 0) {//While the last element of the array is zero:
			if(pp.length == 1) {// If the array only has one element, exit the loop
				break;
			}
			// Create a new array of the same length as the input array
			double [] new_pp = new double[pp.length];

			// Copy the input array into the new array.
			for(int i = 0; i < pp.length;i++) {
				new_pp[i] = pp[i];
			}
			pp = new double[pp.length - 1]; // Create a new array of length one less than the input array.
			// Copy all elements except the last from the new array into the input array.
			for(int i = 0; i < pp.length;i++) {
				pp[i] = new_pp[i];
			}
		}
		return pp;
	}
	/**
	 * This function computes the polynom which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return	
	 * ----------------------------------------
	 * in this function we define a new polynom "ans" and represent the length of the sum of
	 * the length p1 and p2 -1.
	 * after that we go through the coefficients of p1 and p2(started from the last to the first it shpould work also from the fist to the last)
	 *  in every iteration we multiplied the coefficient
	 * in the place [i+j]and adding the  multiplied[i]and[j] 
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = new double [(p1.length+p2.length)-1];

		// Iterate over each coefficient of p1 in reverse order
		for(int i=p1.length-1;i>=0;i--) {

			// Iterate over each coefficient of p2 in reverse order
			for (int j=p2.length-1;j>=0;j--) {

				// Calculate the product of the coefficients
				//and add it to the appropriate index in the answer array
				ans[i+j]+=p1[i]*p2[j];
			}
		}

		return ans;
	}
	/**
	 * This function computes the derivative polynom:.
	 * @param po
	 * @return
	 * 
	 * -----------------------
	 * to compute the derivetive by doing:
	 * f(x)= k*x^n____f'x=n*k*x^n-1
	 * (k should be any number) 
	 * f(x)=k ___ f'(x)= 0
	 */

	public static double[] derivative (double[] po) {
		double[] ans = ZERO;
		int l = po.length;
		if(l>=1) {// If the input array has length 1 or more
			ans = new double[l - 1];
			for (int i = l - 1; i > 0; i--) {//iterates through the po array in from the last to the first

				// Multiply the coefficient by the power
				//to removing the 0 we define ans in [i-1] 
				ans[i - 1] = po[i] * i;    
			}

		}
		return ans;// Return the derivative array
	}

}

