

import org.junit.jupiter.api.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2023, Ariel University,
 *  * Ex2: arrays, static functions and JUnit
 * 
 * This JUnit class represents a JUnit (unit testing) for Ex2 - 
 * It contains few testing functions for the polynum functions as define in Ex2.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex2Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};

	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex2.f(po1, 0);
		double fx1 = Ex2.f(po1, 1);
		double fx2 = Ex2.f(po1, 2);
		assertEquals(fx0, 2, Ex2.EPS);
		assertEquals(fx1, 4, Ex2.EPS);
		assertEquals(fx2, 6, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex2.add(po1, po2);
		double f1x = Ex2.f(po1, x);
		double f2x = Ex2.f(po2, x);
		double f12x = Ex2.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex2.add(po1, po2);//Tests that po1+po2+ (-1*po2) == po1
		double[] minus1 = {-1};
		double[] pp2 = Ex2.mul(po2, minus1);
		double[] p1 = Ex2.add(p12, pp2);
		assertTrue(Ex2.equals(p1, po1));
		double []po5= {-1.0, 2.61, 0.2};

		//Tests that p1+p2+ (-1*p2) == p1
		double [] po6=Ex2.add(P1, P2);//p2+p2
		double [] po7=Ex2.mul(minus1, P2);//p2*-1
		double[] po8=Ex2.add(po6, po7);
		assertTrue(Ex2.equals(po8, P1));
		assertTrue(Ex2.equals(p12, po5));//test that po1+po2=po5

	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		// Tests that po1+po2 == po2+po1
		double[] p12 = Ex2.add(po1, po2);
		double[] p21 = Ex2.add(po2, po1);
		assertTrue(Ex2.equals(p12, p21));

		//Tests that p1+p2 == p2+p1
		double[]p23=Ex2.add(P1, P2);
		double[]p32=Ex2.add(P2, P1);
		assertTrue(Ex2.equals(p23, p32));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		//Tests that po1+0 == po1
		double[] p1 = Ex2.add(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, po1));
		//Tests that p1+0 == p1
		double[] p2 = Ex2.add(p1, Ex2.ZERO);
		assertTrue(Ex2.equals(p2, p1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		//Tests that po1*0 == 0
		double[] p1 = Ex2.mul(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, Ex2.ZERO));
		double[] result=Ex2.mul(po1, po2);
		double[]expected= {-6.0, -4.78, 1.62, 0.4};
		assertArrayEquals(expected,result);

		//Tests that p1*0 == 0
		double[] p2 = Ex2.mul(p1, Ex2.ZERO);
		assertTrue(Ex2.equals(p2, Ex2.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		//Tests that po1*po2 == po2*po1
		double[] p12 = Ex2.mul(po1, po2);
		double[] p21 = Ex2.mul(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
		//Tests that p1*p2 == p2*p1	
		double[] xx = Ex2.mul(po1, po2);
		double[] yy = Ex2.mul(po2, po1);
		assertTrue(Ex2.equals(p12, p21));

	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex2.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex2.f(po1, x);
			double f2x = Ex2.f(po2, x);
			double f12x = Ex2.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex2.EPS);
			double[] p23 = Ex2.mul(P1, P2);
			//Tests that p1(x) * p2(x) = (p1*p2)(x),
			double f11x = Ex2.f(P1, x);
			double f22x = Ex2.f(P2, x);
			double f121x = Ex2.f(p23, x);
			assertEquals(f121x, f11x*f22x, Ex2.EPS);
			//test the result of  p1*p2
			double[]expexted222= {0.2, 0.0, 2.3, 0.1, 9.0, -0.7, 8.9, -3.0, 0.0};
			assertArrayEquals(expexted222,p23);

		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex2.derivative(p); // 6x + 2
		double[] dp2 = Ex2.derivative(dp1); // 6
		double[] dp3 = Ex2.derivative(dp2); // 0
		double[] dp4 = Ex2.derivative(dp3); // 0

		assertTrue(Ex2.equals(dp1, pt));//check by use the equals function if the derivative(p)
		//is equal to the polynom pt.

		assertTrue(Ex2.equals(Ex2.ZERO, dp3));//check if dp3 equal to zero.
		assertEquals(dp1[1], pt[1],Ex2.EPS);//check if the coffecient in the same index is the same
		assertEquals(1,dp2.length);
		assertTrue(Ex2.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";

		double[]result= Ex2.getPolynomFromString(sp2);
		assertTrue(Ex2.equals(p, result));

		String sp = Ex2.poly(p);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex2.ZERO, {1+Ex2.EPS/2}, {1,2}};
		double[][] xx = {{-2*Ex2.EPS}, {1+Ex2.EPS*1.2}, {1,2,Ex2.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex2.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex2.equals(d1[i], xx[i]));

		}
		// i added another cases to test
		double [] pp1= {1,2,3};
		double [] pp2= {1,2,3};
		double [] pp3= {5.6,7};
		double [] pp4= {0,0,2,0};
		double [] pp5= {0,0,2,0,0};

		boolean expected= Ex2.equals(po1, po2);
		boolean expected1= Ex2.equals(pp1, pp2);
		boolean expected2= Ex2.equals(pp2, pp3);
		boolean expected3= Ex2.equals(pp4, pp5);
		assertFalse(expected);
		assertTrue(expected1);
		assertFalse(expected2);
		assertTrue(expected3);

	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue() {
		double x1=0, x2=4;
		double rs1 = Ex2.sameValue(po1,po2, x1, x2, Ex2.EPS);
		double rs2 = Ex2.sameValue(po2,po1, x1, x2, Ex2.EPS);
		assertEquals(rs1,rs2,Ex2.EPS);

	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=0, x2=4;
		double a1 = Ex2.area(po1, po2, x1, x2, 1000);
		double a2 = Ex2.area(po2, po1, x1, x2, 1000);
		assertEquals(a1,a2,Ex2.EPS);
		// test another case
		double[] xx= {0,-1};
		double[] xx1= {0,0,1};
		double x_1=0,x_2=-1;
		double result =Ex2.area(xx, xx1, x_2, x_1, 1000);
		assertEquals(0.16666650000000005,result);
	}

	@Test
	/*test the String Poly function
	 * in this test we test two conditions
	 * 1. test p1//
	 * 2.test po2
	 */
	public void testPoly() {
		String result1=Ex2.poly(P1);
		assertEquals("-x^3 + 3.0x^2 + 2.0",result1);
		String result2= Ex2.poly(po2);
		assertEquals("0.2x^2 + 0.61x - 3.0",result2);
		String result3=Ex2.poly(po1);//2.0x+2.0

		double[] arr= {2,2,4,3};
		String result4= Ex2.poly(arr);
		String expected="3.0x^3 + 4.0x^2 + 2.0x + 2.0";
		assertEquals(expected,result4);

	}
	/*
	 * this test of the poly function
	 * i test the case with thre points and
	 * the case with two points
	 */
	@Test
	public void testPolynomFromPoints() {
		double[] xx = {1, 2, 3};
		double[] yy = {1, 4, 9};
		double[] result = Ex2.PolynomFromPoints(xx, yy);
		double[] expected = {0.0, 0.0, 1.0};

		assertArrayEquals(expected, result, Ex2.EPS);

		double[] xx1 = {1,2};
		double[] yy1 = {1,3};
		double[] result1 = Ex2.PolynomFromPoints(yy1, xx1);
		double[] expected1 = {0.5, 0.5};
		assertArrayEquals(expected1, result1, Ex2.EPS);

		double[] xx2 = {1, 2, 3};
		double[] yy2 = {1, 2, 4};
		double[] result2 = Ex2.PolynomFromPoints(xx, yy);
		double[] expected2 = {1.0, -0.5, 0.5};

		assertArrayEquals(expected, result, Ex2.EPS);
	}
	@Test
	public void testlength() {
		Double result= Ex2.length(po1, 0, 6, 3)	;
		double expected=3*Math.sqrt(20)	;
		assertEquals(expected,result);

	}
}


