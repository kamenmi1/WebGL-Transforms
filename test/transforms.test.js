"use strict";

describe("ZeroArray", () => {

	it("works with parameter", () => {
		let x = new ZeroArray(2);
		expect(x.length).toBe(2);
		expect(x[0].length).toBe(2);
		expect(x[1].length).toBe(2);
	});

	it("returns 4-dimensional array when called w/o parameter", () => {
		let x = new ZeroArray();
		expect(x.length).toBe(4);
		expect(x[2].length).toBe(4);
		expect(x[3].length).toBe(4);
	});

	it("throws TypeError when parameter is not a number", () => {
		expect(() => new ZeroArray("a")).toThrowError(TypeError);
	});
});

describe("Vec1D", () => {

	it("constructor can make a copy of another Vec1D", () => {
		let v = new Vec1D(new Vec1D(5));
		expect(v.x).toBe(5);
	});

	it("constructor w/o parameter sets 0", () => {
		let v = new Vec1D();
		expect(v.x).toBe(0);
	});

	it("constructor w/ parameter works", () => {
		let v = new Vec1D(1);
		expect(v.x).toBe(1);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Vec1D("a")).toThrowError(TypeError);
	});

	it("add prototype works with Vec1D as parameter", () => {
		let v = new Vec1D(1);
		let w = v.add(new Vec1D(3));
		expect(w.x).toBe(4);
	});

	it("add prototype throws TypeError when parameter is not a Vec1D", () => {
		expect(() => new Vec1D().add("a")).toThrowError(TypeError);
		expect(() => new Vec1D().add()).toThrowError(TypeError);
	});

	it("mul prototype works with number as parameter", () => {
		let v = new Vec1D(2);
		let w = v.mul(3);
		expect(w.x).toBe(6);
	});

	it("mul prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Vec1D().mul("a")).toThrowError(TypeError);
		expect(() => new Vec1D().mul()).toThrowError(TypeError);
	});

	it("c function returns itself", () => {
		expect(new Vec1D().c() instanceof Vec1D).toBe(true);
	});
});

describe("Vec2D", () => {

	it("constructor can make a copy of another Vec2D", () => {
		let v = new Vec2D(new Vec2D(5, 6));
		expect(v.x).toBe(5);
		expect(v.y).toBe(6);
	});

	it("constructor w/o parameters sets 0", () => {
		let v = new Vec2D();
		expect(v.x).toBe(0);
		expect(v.y).toBe(0);
	});

	it("constructor w/ parameters works", () => {
		let v = new Vec2D(-2, -3);
		expect(v.x).toBe(-2);
		expect(v.y).toBe(-3);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Vec2D("a")).toThrowError(TypeError);
	});

	it("add prototype works with Vec2D as parameter", () => {
		let v = new Vec2D(1, 2);
		let w = v.add(new Vec2D(7, 4));
		expect(w.x).toBe(8);
		expect(w.y).toBe(6);
	});

	it("add prototype throws TypeError when parameter is not a Vec2D", () => {
		expect(() => new Vec2D("a")).toThrowError(TypeError);
		expect(() => new Vec2D(1, true)).toThrowError(TypeError);
		expect(() => new Vec2D(new Vec2D(), "a")).not.toThrowError(TypeError);
	});

	it("mul prototype works with Vec2D as parameter", () => {
		let v = new Vec2D(2, -4);
		let w = v.mul(new Vec2D(5, 7));
		expect(w.x).toBe(10);
		expect(w.y).toBe(-28);
	});

	it("mul prototype works with number as parameter", () => {
		let v = new Vec2D(2, -4);
		let w = v.mul(5);
		expect(w.x).toBe(10);
		expect(w.y).toBe(-20);
	});

	it("mul prototype throws TypeError when parameter is not a number or Vec2D", () => {
		expect(() => new Vec2D().mul("a")).toThrowError(TypeError);
		expect(() => new Vec2D().mul()).toThrowError(TypeError);
	});

	it("dot prototype works with Vec2D as parameter", () => {
		let v = new Vec2D(5, -1);
		let w = v.dot(new Vec2D(2, 6));
		expect(w).toBe(4);
	});

	it("dot prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Vec2D().dot("a")).toThrowError(TypeError);
		expect(() => new Vec2D().dot()).toThrowError(TypeError);
	});

	it("length prototype works as expected", () => {
		let v = new Vec2D(4, 3);
		let l = v.length();
		expect(l).toBe(5);
	});

	it("normalized prototype works as expected", () => {
		let v = new Vec2D(4, 3);
		let w = v.normalized();
		expect(w.x).toBe(0.8);
		expect(w.y).toBe(0.6);

		let t = new Vec3D(0, 0, 0).normalized();
		expect(t.x).toBe(0);
		expect(t.y).toBe(0);
	});

	it("c function returns itself", () => {
		expect(new Vec2D().c() instanceof Vec2D).toBe(true);
	});
});

describe("Vec3D", () => {

	it("constructor can make a copy of another Vec3D", () => {
		let v = new Vec3D(new Vec3D(5, 6, 7));
		expect(v.x).toBe(5);
		expect(v.y).toBe(6);
		expect(v.z).toBe(7);
	});

	it("constructor w/o parameters sets 0", () => {
		let v = new Vec3D();
		expect(v.x).toBe(0);
		expect(v.y).toBe(0);
		expect(v.z).toBe(0);
	});

	it("constructor w/ parameters works", () => {
		let v = new Vec3D(-2, -3, -9);
		expect(v.x).toBe(-2);
		expect(v.y).toBe(-3);
		expect(v.z).toBe(-9);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Vec3D("a")).toThrowError(TypeError);
		expect(() => new Vec3D(1, 1, true)).toThrowError(TypeError);
		expect(() => new Vec3D(new Vec3D(), "a", 2)).not.toThrowError(TypeError);
	});

	it("add prototype works with Vec3D as parameter", () => {
		let v = new Vec3D(1, 2, 5);
		let w = v.add(new Vec3D(7, 4, -5));
		expect(w.x).toBe(8);
		expect(w.y).toBe(6);
		expect(w.z).toBe(0);
	});

	it("add prototype throws TypeError when parameter is not a Vec3D", () => {
		expect(() => new Vec3D().add("a")).toThrowError(TypeError);
		expect(() => new Vec3D().add()).toThrowError(TypeError);
	});

	it("mul prototype works with Vec3D as parameter", () => {
		let v = new Vec3D(2, -4, 3.5);
		let w = v.mul(new Vec3D(5, 7, -1));
		expect(w.x).toBe(10);
		expect(w.y).toBe(-28);
		expect(w.z).toBe(-3.5);
	});

	it("mul prototype works with number as parameter", () => {
		let v = new Vec3D(2, -4, 1);
		let w = v.mul(-2);
		expect(w.x).toBe(-4);
		expect(w.y).toBe(8);
		expect(w.z).toBe(-2);
	});

	it("mul prototype works with Mat3 as parameter", () => {
		let m = new Mat3(new Vec3D(1, 2, 3), new Vec3D(-2, 0, 3), new Vec3D(1.5, 8, -1));
		let v = new Vec3D(3, -2, 3);
		let w = v.mul(m);
		expect(w.x).toBe(11.5);
		expect(w.y).toBe(30);
		expect(w.z).toBe(0);
	});
/*
	it("mul prototype works with Quat as parameter", () => {
		// TODO
	});
*/
	it("mul prototype throws TypeError when parameter is not a number or Vec3D or Mat3 or Quat", () => {
		expect(() => new Vec3D().mul("a")).toThrowError(TypeError);
		expect(() => new Vec3D().mul()).toThrowError(TypeError);
	});

	it("dot prototype works with Vec3D as parameter", () => {
		let v = new Vec3D(5, -1, 2);
		let d = v.dot(new Vec3D(7, 2, 6));
		expect(d).toBe(45);
	});

	it("dot prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Vec3D().dot("a")).toThrowError(TypeError);
		expect(() => new Vec3D().dot()).toThrowError(TypeError);
	});

	it("cross prototype works with Vec3D as parameter", () => {
		let v = new Vec3D(5, -1, 2);
		let w = v.cross(new Vec3D(7, 2, 6));
		expect(w.x).toBe(-10);
		expect(w.y).toBe(-16);
		expect(w.z).toBe(17);
	});

	it("cross prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Vec3D().cross("a")).toThrowError(TypeError);
		expect(() => new Vec3D().cross()).toThrowError(TypeError);
	});

	it("length prototype works as expected", () => {
		let v = new Vec3D(4, 0, 3);
		let l = v.length();
		expect(l).toBe(5);
	});

	it("normalized prototype works as expected", () => {
		let v = new Vec3D(4, 0, 3);
		let w = v.normalized();
		expect(w.x).toBe(0.8);
		expect(w.y).toBe(0);
		expect(w.z).toBe(0.6);

		let t = new Vec3D(0, 0, 0).normalized();
		expect(t.x).toBe(0);
		expect(t.y).toBe(0);
		expect(t.z).toBe(0);
	});

	it("c function returns itself", () => {
		expect(new Vec3D().c() instanceof Vec3D).toBe(true);
	});
});


describe("Point3D", () => {

	it("constructor can make a copy of another Point3D", () => {
		let v = new Point3D(new Point3D(5, 6, 7, 8));
		expect(v.x).toBe(5);
		expect(v.y).toBe(6);
		expect(v.z).toBe(7);
		expect(v.w).toBe(8);
	});

	it("constructor can make a copy of another Vec3D", () => {
		let v = new Point3D(new Vec3D(5, 7, 8));
		expect(v.x).toBe(5);
		expect(v.y).toBe(7);
		expect(v.z).toBe(8);
		expect(v.w).toBe(1);
	});

	it("constructor w/o parameters sets 0", () => {
		let v = new Point3D();
		expect(v.x).toBe(0);
		expect(v.y).toBe(0);
		expect(v.z).toBe(0);
		expect(v.w).toBe(1);
	});

	it("constructor w/ parameters works", () => {
		let v = new Point3D(-2, -3, -9, 3);
		expect(v.x).toBe(-2);
		expect(v.y).toBe(-3);
		expect(v.z).toBe(-9);
		expect(v.w).toBe(3);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Point3D("a")).toThrowError(TypeError);
		expect(() => new Point3D(1, 1, true)).toThrowError(TypeError);
		expect(() => new Point3D(new Point3D(), "a", 2)).not.toThrowError(TypeError);
	});

	it("add prototype works with Point3D as parameter", () => {
		let v = new Point3D(-4, 1, 2, 5);
		let w = v.add(new Point3D(1.3, 7, 4, -5));
		expect(w.x).toBe(-2.7);
		expect(w.y).toBe(8);
		expect(w.z).toBe(6);
		expect(w.w).toBe(0);
	});

	it("add prototype throws TypeError when parameter is not a Point3D", () => {
		expect(() => new Point3D().add("a")).toThrowError(TypeError);
		expect(() => new Point3D().add()).toThrowError(TypeError);
	});

	it("mul prototype works with number as parameter", () => {
		let v = new Point3D(2, -4, 1, 2.5);
		let w = v.mul(-2);
		expect(w.x).toBe(-4);
		expect(w.y).toBe(8);
		expect(w.z).toBe(-2);
		expect(w.w).toBe(-5);
	});

	it("mul prototype works with Mat4 as parameter", () => {
		let m = new Mat4(new Point3D(1, 2, 3, 4), new Point3D(-1, -9, 5, 2), new Point3D(3, 0, 1, -2), new Point3D(5, 1, -3, 1));
		let v = new Point3D(1, -4, 1.5, 2);
		let w = v.mul(m);
		expect(w.x).toBe(19.5);
		expect(w.y).toBe(40);
		expect(w.z).toBe(-21.5);
		expect(w.w).toBe(-5);
	});
/*
	it("mul prototype works with Quat as parameter", () => {
		// TODO
	});
*/
	it("mul prototype throws TypeError when parameter is not a number or Mat4 or Quat", () => {
		expect(() => new Point3D().mul("a")).toThrowError(TypeError);
		expect(() => new Point3D().mul()).toThrowError(TypeError);
	});

	it("sub prototype works with Point3D as parameter", () => {
		let v = new Point3D(5, -1, 2, 3);
		let w = v.sub(new Point3D(-2, 7, 2, 6));
		expect(w.x).toBe(7);
		expect(w.y).toBe(-8);
		expect(w.z).toBe(0);
		expect(w.w).toBe(-3);
	});

	it("sub prototype throws TypeError when parameter is not a Point3D", () => {
		expect(() => new Point3D().dot("a")).toThrowError(TypeError);
		expect(() => new Point3D().dot()).toThrowError(TypeError);
	});

	it("dehomog prototype returns correct Vec3D", () => {
		let v1 = new Point3D(3, 5, -1, 2).dehomog();
		expect(v1 instanceof Vec3D).toBe(true);
		expect(v1.x).toBe(1.5);
		expect(v1.y).toBe(2.5);
		expect(v1.z).toBe(-0.5);
		let v2 = new Point3D(3, 5, -1, 0).dehomog();
		expect(v2 instanceof Vec3D).toBe(true);
		expect(v2.x).toBe(0);
		expect(v2.y).toBe(0);
		expect(v2.z).toBe(0);
	});

	it("ignoreW prototype returns correct Vec3D", () => {
		let v = new Point3D(4, 0, 3, 1).ignoreW();
		expect(v instanceof Vec3D).toBe(true);
		expect(v.x).toBe(4);
		expect(v.y).toBe(0);
		expect(v.z).toBe(3);
	});

	it("c function returns itself", () => {
		expect(new Point3D().c() instanceof Point3D).toBe(true);
	});
});


describe("Quat", () => {

	it("constructor can make a copy of another Quat", () => {
		let q = new Quat(new Quat(5, 6, 7, 8));
		expect(q.r).toBe(5);
		expect(q.i).toBe(6);
		expect(q.j).toBe(7);
		expect(q.k).toBe(8);
	});

	it("constructor w/o parameters sets 0", () => {
		let q = new Quat();
		expect(q.r).toBe(0);
		expect(q.i).toBe(0);
		expect(q.j).toBe(0);
		expect(q.k).toBe(0);
	});

	it("constructor w/ parameters works", () => {
		let q = new Quat(-2, -3, -9, 3);
		expect(q.r).toBe(-2);
		expect(q.i).toBe(-3);
		expect(q.j).toBe(-9);
		expect(q.k).toBe(3);
	});

	it("constructor throws TypeError when parameter is not a number or Quat", () => {
		expect(() => new Quat("a")).toThrowError(TypeError);
		expect(() => new Quat(1, 1, true)).toThrowError(TypeError);
		expect(() => new Quat(new Quat(), "a", 2)).not.toThrowError(TypeError);
	});

	it("add prototype works with Quat as parameter", () => {
		let v = new Quat(-4, 1, 2, 5);
		let q = v.add(new Quat(1.3, 7, 4, -5));
		expect(q.r).toBe(-2.7);
		expect(q.i).toBe(8);
		expect(q.j).toBe(6);
		expect(q.k).toBe(0);
	});

	it("add prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().add("a")).toThrowError(TypeError);
		expect(() => new Quat().add()).toThrowError(TypeError);
	});

	it("sub prototype works with Quat as parameter", () => {
		let v = new Quat(-4, 1, 2, 5);
		let q = v.sub(new Quat(1.3, 7, 4, -5));
		expect(q.r).toBe(-5.3);
		expect(q.i).toBe(-6);
		expect(q.j).toBe(-2);
		expect(q.k).toBe(10);
	});

	it("sub prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().sub("a")).toThrowError(TypeError);
		expect(() => new Quat().sub()).toThrowError(TypeError);
	});

	it("mulR prototype works with Quat as parameter", () => {
		let v = new Quat(1, 3, -2, 7);
		let q = v.mulR(new Quat(-2, 0, 6, 3));
		expect(q.r).toBe(-11);
		expect(q.i).toBe(-54);
		expect(q.j).toBe(1);
		expect(q.k).toBe(7);
	});

	it("mulR prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().mulR("a")).toThrowError(TypeError);
		expect(() => new Quat().mulR()).toThrowError(TypeError);
	});

	it("mulL prototype works with Quat as parameter", () => {
		let v = new Quat(1, 3, -2, 7);
		let q = v.mulL(new Quat(-2, 0, 6, 3));
		expect(q.r).toBe(-11);
		expect(q.i).toBe(42);
		expect(q.j).toBe(19);
		expect(q.k).toBe(-29);
	});

	it("mulL prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().mulL("a")).toThrowError(TypeError);
		expect(() => new Quat().mulL()).toThrowError(TypeError);
	});

	it("mul prototype works with Quat as parameter", () => {
		let v = new Quat(1, 3, -2, 7);
		let q = v.mul(new Quat(-2, 0, 6, 3));
		expect(q.r).toBe(-11);
		expect(q.i).toBe(-54);
		expect(q.j).toBe(1);
		expect(q.k).toBe(7);
	});

	it("mul prototype works with number as parameter", () => {
		let v = new Quat(1, 3, -2, 7);
		let q = v.mul(3);
		expect(q.r).toBe(3);
		expect(q.i).toBe(9);
		expect(q.j).toBe(-6);
		expect(q.k).toBe(21);
	});

	it("mul prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().mul("a")).toThrowError(TypeError);
		expect(() => new Quat().mul()).toThrowError(TypeError);
	});

	it("norma prototype works as expected", () => {
		let n = new Quat(1, 6, -4, 7).norma();
		expect(n).toBeCloseTo(10.099504938362077, 5);
	});

	it("inv prototype works as expected", () => {
		let q1 = new Quat(1, 6, -4, 7).inv();
		expect(q1.r).toBeCloseTo(0.009803921568627453, 5);
		expect(q1.i).toBeCloseTo(-0.05882352941176471, 5);
		expect(q1.j).toBeCloseTo(0.03921568627450981, 5);
		expect(q1.k).toBeCloseTo(-0.06862745098039216, 5);

		let q2 = new Quat().inv();
		expect(q2.r).toBe(0);
		expect(q2.i).toBe(0);
		expect(q2.j).toBe(0);
		expect(q2.k).toBe(0);
	});

	it("log prototype works as expected", () => {
		let q1 = new Quat(1, 3, -2, 7).log();
		expect(q1.r).toBeCloseTo(2.0715673631957663, 5);
		expect(q1.i).toBeCloseTo(0.5503445589964804, 5);
		expect(q1.j).toBeCloseTo(-0.3668963726643203, 5);
		expect(q1.k).toBeCloseTo(1.284137304325121, 5);

		let q2 = new Quat(-2, 0, 0, 0).log();
		expect(q2.r).toBeCloseTo(0.6931471805599453, 5);
		expect(q2.i).toBe(1);
		expect(q2.j).toBe(0);
		expect(q2.k).toBe(0);

		let q3 = new Quat(2, 0, 0, 0).log();
		expect(q3.r).toBeCloseTo(0.6931471805599453, 5);
		expect(q3.i).toBe(0);
		expect(q3.j).toBe(0);
		expect(q3.k).toBe(0);

		let q4 = new Quat().log();
		expect(q4.r).toBe(0);
		expect(q4.i).toBe(0);
		expect(q4.j).toBe(0);
		expect(q4.k).toBe(0);
	});

	it("exp prototype works as expected", () => {
		let q1 = new Quat(1, 3, -2, 7).exp();
		expect(q1.r).toBeCloseTo(-0.05443332579873212, 5);
		expect(q1.i).toBeCloseTo(1.0354587420929653, 5);
		expect(q1.j).toBeCloseTo(-0.6903058280619768, 5);
		expect(q1.k).toBeCloseTo(2.416070398216919, 5);

		let q2 = new Quat(-2, 0, 0, 0).exp();
		expect(q2.r).toBeCloseTo(0.1353352832366127, 5);
		expect(q2.i).toBe(0);
		expect(q2.j).toBe(0);
		expect(q2.k).toBe(0);
	});

	it("neg prototype works as expected", () => {
		let q1 = new Quat(1, 3, -2, 7).neg();
		expect(q1.r).toBe(-1);
		expect(q1.i).toBe(-3);
		expect(q1.j).toBe(2);
		expect(q1.k).toBe(-7);
	});

	it("dot prototype works as expected", () => {
		let q1 = new Quat(1, 3, -2, 7);
		let q2 = new Quat(2, -6, 1.5, 1);
		let dot = q1.dot(q2);
		expect(dot).toBe(-12);
	});

	it("dot prototype throws TypeError when parameter is not a Quat", () => {
		expect(() => new Quat().dot("a")).toThrowError(TypeError);
		expect(() => new Quat().dot()).toThrowError(TypeError);
	});

	it("renorm prototype works as expected", () => {
		let q1 = new Quat(1, 3, -2, 7).renorm();
		expect(q1.r).toBeCloseTo(0.1259881576697424, 5);
		expect(q1.i).toBeCloseTo(0.3779644730092272, 5);
		expect(q1.j).toBeCloseTo(-0.2519763153394848, 5);
		expect(q1.k).toBeCloseTo(0.8819171036881969, 5);
	});

	it("toRotationMatrix prototype works as expected", () => {
		let m = new Quat(1, 3, -2, 7).toRotationMatrix().mat;
		expect(m[0][0]).toBe(-105);
		expect(m[0][1]).toBe(2);
		expect(m[0][2]).toBe(46);
		expect(m[0][3]).toBe(0);
		expect(m[1][0]).toBe(-26);
		expect(m[1][1]).toBe(-115);
		expect(m[1][2]).toBe(-22);
		expect(m[1][3]).toBe(0);
		expect(m[2][0]).toBe(38);
		expect(m[2][1]).toBe(-34);
		expect(m[2][2]).toBe(-25);
		expect(m[2][3]).toBe(0);
		expect(m[3][0]).toBe(0);
		expect(m[3][1]).toBe(0);
		expect(m[3][2]).toBe(0);
		expect(m[3][3]).toBe(1);
	});

	it("toEulerAngle prototype works as expected", () => {
		let v = new Quat(0.4, 3, -2, 7).toEulerAngle();
		expect(v.x).toBeCloseTo(2.318558961454817, 5);
		expect(v.y).toBeCloseTo(0.3810003810005715, 5);
		expect(v.z).toBeCloseTo(-0.254000254000381, 5);
		expect(v.w).toBeCloseTo(0.8890008890013334, 5);
	});

	it("lerp prototype works as expected", () => {
		let q1 = new Quat(0.4, 3, -2, 7);
		let q2 = new Quat(4, -2, 3, 7);

		let q3 = q1.lerp(q2, 2);
		expect(q3.r).toBe(q2.r);
		expect(q3.i).toBe(q2.i);
		expect(q3.j).toBe(q2.j);
		expect(q3.k).toBe(q2.k);

		let q4 = q1.lerp(q2, 0);
		expect(q4.r).toBe(q1.r);
		expect(q4.i).toBe(q1.i);
		expect(q4.j).toBe(q1.j);
		expect(q4.k).toBe(q1.k);

		let q5 = q1.lerp(q2, 0.4);
		let q6 = q2.lerp(q1, 0.6);
		expect(q5.r).toBe(q6.r);
		expect(q5.i).toBe(q6.i);
		expect(q5.j).toBe(q6.j);
		expect(q5.k).toBe(q6.k);
	});

	it("lerp prototype throws TypeError when parameters are not Quat and number", () => {
		expect(() => new Quat().lerp("a", "a")).toThrowError(TypeError);
		expect(() => new Quat().lerp()).toThrowError(TypeError);
	});

	it("slerp prototype works as expected", () => {
		let q1 = new Quat(5, 3, -2, 7);

		let q2 = new Quat(1, 2, -5, -3);// dot 0
		let q3 = q1.slerp(q2, 0.5);
		expect(q3.r).toBeCloseTo(0.49227672520391125, 5);
		expect(q3.i).toBeCloseTo(0.45388481990259016, 5);
		expect(q3.j).toBeCloseTo(-0.7177581257880786, 5);
		expect(q3.k).toBeCloseTo(0.19098552026185375, 5);

		let q4 = new Quat(1, 2, -3, -1);// dot 10
		let q5 = q1.slerp(q4, 0);
		expect(q5.r).toBe(5);
		expect(q5.i).toBe(3);
		expect(q5.j).toBe(-2);
		expect(q5.k).toBe(7);

		let q6 = new Quat(1, 2, -3, -3);// dot -4
		let q7 = q1.slerp(q6, 0.4);
		expect(q7.r).toBeCloseTo(0.5517368198016962, 5);
		expect(q7.i).toBeCloseTo(0.547358844043825, 5);
		expect(q7.j).toBeCloseTo(-0.6224258390801777, 5);
		expect(q7.k).toBeCloseTo(0.09257889806783776, 5);

		let q8 = new Quat(1, 2, -5, -3.1);// dot -0.7
		let q9 = q1.slerp(q8, 2);
		expect(q9.r).toBe(q1.r);
		expect(q9.i).toBe(q1.i);
		expect(q9.j).toBe(q1.j);
		expect(q9.k).toBe(q1.k);

		let q10 = q1.slerp(q8, -1);
		expect(q10.r).toBe(q8.r);
		expect(q10.i).toBe(q8.i);
		expect(q10.j).toBe(q8.j);
		expect(q10.k).toBe(q8.k);
	});

	it("slerp prototype throws TypeError when parameters are not Quat and number", () => {
		expect(() => new Quat().slerp("a", "a")).toThrowError(TypeError);
		expect(() => new Quat().slerp()).toThrowError(TypeError);
	});

	it("squad prototype works as expected", () => {
		let a = new Quat(5, 3, -2, 7);
		let b = new Quat(1, 2, -5, -3.1);
		let c = new Quat(-2, 0, 6, 3);
		let d = new Quat(-1, 0, 1, 2);
		let e = a.squad(b, c, d, 0.6);

		expect(e.r).toBeCloseTo(0.24857015719081804, 5);
		expect(e.i).toBeCloseTo(0.6269395906460726, 5);
		expect(e.j).toBeCloseTo(0.14129984737341578, 5);
		expect(e.k).toBeCloseTo(0.7247026837034042, 5);
	});

	it("squad prototype throws TypeError when parameters are not Quat-s and number", () => {
		expect(() => new Quat().squad("a", "a", "a", "a")).toThrowError(TypeError);
		expect(() => new Quat().squad()).toThrowError(TypeError);
	});

	it("quadrangle prototype works as expected", () => {
		let a = new Quat(5, 3, -2, 7);
		let b = new Quat(1, 2, -5, -3.1);
		let c = new Quat(-2, 0, 6, 3);
		let e = a.quadrangle(b, c);

		expect(e.r).toBeCloseTo(1.1738922180666038, 5);
		expect(e.i).toBeCloseTo(-0.11515063497536307, 5);
		expect(e.j).toBeCloseTo(0.11758452348157551, 5);
		expect(e.k).toBeCloseTo(0.012650862133917222, 5);
	});

	it("quadrangle prototype throws TypeError when parameters are not Quat-s and number", () => {
		expect(() => new Quat().quadrangle("a", "a")).toThrowError(TypeError);
		expect(() => new Quat().quadrangle()).toThrowError(TypeError);
	});

	it("squad2 prototype works as expected", () => {
		let a = new Quat(5, 3, -2, 7);
		let b = new Quat(1, 2, -5, -3.1);
		let c = new Quat(-2, 0, 6, 3);
		let d = new Quat(-1, 0, 1, 2);
		let e = a.squad2(b, c, d, 0.6);

		expect(e.r).toBeCloseTo(0.73357337319017, 5);
		expect(e.i).toBeCloseTo(0.09203478946363089, 5);
		expect(e.j).toBeCloseTo(0.3660624883267518, 5);
		expect(e.k).toBeCloseTo(0.5651530397288929, 5);
	});

	it("squad2 prototype throws TypeError when parameters are not Quat-s and number", () => {
		expect(() => new Quat().squad2("a", "a", "a", "a")).toThrowError(TypeError);
		expect(() => new Quat().squad2()).toThrowError(TypeError);
	});

	it("c function returns itself", () => {
		expect(new Quat().c() instanceof Quat).toBe(true);
	});
});


describe("Quat2", () => {

	it("fromRotationMatrix works as expected", () => {
		let m = new Mat4(new Point3D(1, 2, 3, 4), new Point3D(-1, -9, 5, 2), new Point3D(3, 0, 1, -2), new Point3D(5, 1, -3, 1));
		let q = Quat2.fromRotationMatrix(m);

		expect(q.r).toBeCloseTo(1.5811388300841898, 5);
		expect(q.i).toBeCloseTo(0.15811388300841897, 5);
		expect(q.j).toBeCloseTo(0.9486832980505138, 5);
		expect(q.k).toBeCloseTo(0.7905694150420948, 5);
	});

	it("fromRotationMatrix throws TypeError when parameter is not Mat4", () => {
		expect(() => Quat2.fromRotationMatrix("a")).toThrowError(TypeError);
		expect(() => Quat2.fromRotationMatrix()).toThrowError(TypeError);
	});

	it("fromEulerAngle works as expected", () => {
		let q = Quat2.fromEulerAngle(45, 0.5, 1, 0);

		expect(q.r).toBeCloseTo(-0.8733046400935156, 5);
		expect(q.i).toBeCloseTo(-0.24358725623025476, 5);
		expect(q.j).toBeCloseTo(-0.4871745124605095, 5);
		expect(q.k).toBeCloseTo(-0, 5);
	});

	it("fromEulerAngle throws TypeError when parameters are not numbers", () => {
		expect(() => Quat2.fromEulerAngle("a", "a", "a", "a")).toThrowError(TypeError);
		expect(() => Quat2.fromEulerAngle()).toThrowError(TypeError);
	});

	it("fromEulerAngles works as expected", () => {
		let q = Quat2.fromEulerAngles(45, 60, 125);

		expect(q.r).toBeCloseTo(-0.28417852247107306, 5);
		expect(q.i).toBeCloseTo(0.21006619016961586, 5);
		expect(q.j).toBeCloseTo(0.8402582556996262, 5);
		expect(q.k).toBeCloseTo(-0.4111943905775006, 5);
	});

	it("fromEulerAngles throws TypeError when parameters are not numbers", () => {
		expect(() => Quat2.fromEulerAngles("a", "a", "a")).toThrowError(TypeError);
		expect(() => Quat2.fromEulerAngles()).toThrowError(TypeError);
	});
});


describe("Mat3", () => {
	let v1 = new Vec3D(4, -2, 3);
	let v2 = new Vec3D(1, 2, -3);
	let v3 = new Vec3D(0, -1, 7);
	let m1 = new Mat3(v1, v2, v3);
	let v4 = new Vec3D(1, -5, 1);
	let v5 = new Vec3D(-2, 2, 4);
	let v6 = new Vec3D(3, -5, 1);
	let m2 = new Mat3(v4, v5, v6);

	it("constructor works correctly with Vec3D as parameters", () => {
		expect(m1.mat[0][0]).toBe(4);
		expect(m1.mat[0][1]).toBe(-2);
		expect(m1.mat[0][2]).toBe(3);
		expect(m1.mat[1][0]).toBe(1);
		expect(m1.mat[1][1]).toBe(2);
		expect(m1.mat[1][2]).toBe(-3);
		expect(m1.mat[2][0]).toBe(0);
		expect(m1.mat[2][1]).toBe(-1);
		expect(m1.mat[2][2]).toBe(7);
	});

	it("constructor works correctly with Mat3 as parameter", () => {
		let m3 = new Mat3(m1);

		expect(m3.mat[0][0]).toBe(4);
		expect(m3.mat[0][1]).toBe(-2);
		expect(m3.mat[0][2]).toBe(3);
		expect(m3.mat[1][0]).toBe(1);
		expect(m3.mat[1][1]).toBe(2);
		expect(m3.mat[1][2]).toBe(-3);
		expect(m3.mat[2][0]).toBe(0);
		expect(m3.mat[2][1]).toBe(-1);
		expect(m3.mat[2][2]).toBe(7);
	});

	it("constructor works correctly with Mat4 as parameter", () => {
		let m3 = new Mat4Identity();
		let m4 = new Mat3(m3);

		expect(m4.mat[0][0]).toBe(1);
		expect(m4.mat[0][1]).toBe(0);
		expect(m4.mat[0][2]).toBe(0);
		expect(m4.mat[1][0]).toBe(0);
		expect(m4.mat[1][1]).toBe(1);
		expect(m4.mat[1][2]).toBe(0);
		expect(m4.mat[2][0]).toBe(0);
		expect(m4.mat[2][1]).toBe(0);
		expect(m4.mat[2][2]).toBe(1);
	});

	it("constructor works correctly w/o parameters", () => {
		let m = new Mat3();

		expect(m.mat[0][0]).toBe(0);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(0);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(0);
	});

	it("add prototype works with Mat3 as parameter", () => {
		let m3 = m1.add(m2);

		expect(m3.mat[0][0]).toBe(5);
		expect(m3.mat[0][1]).toBe(-7);
		expect(m3.mat[0][2]).toBe(4);
		expect(m3.mat[1][0]).toBe(-1);
		expect(m3.mat[1][1]).toBe(4);
		expect(m3.mat[1][2]).toBe(1);
		expect(m3.mat[2][0]).toBe(3);
		expect(m3.mat[2][1]).toBe(-6);
		expect(m3.mat[2][2]).toBe(8);
	});

	it("add prototype throws TypeError when parameter is not a Mat3", () => {
		expect(() => new Mat3().add("a")).toThrowError(TypeError);
		expect(() => new Mat3().add()).toThrowError(TypeError);
	});

	it("mul prototype works with Mat3 as parameter", () => {
		let m3 = m1.mul(m2);

		expect(m3.mat[0][0]).toBe(17);
		expect(m3.mat[0][1]).toBe(-39);
		expect(m3.mat[0][2]).toBe(-1);
		expect(m3.mat[1][0]).toBe(-12);
		expect(m3.mat[1][1]).toBe(14);
		expect(m3.mat[1][2]).toBe(6);
		expect(m3.mat[2][0]).toBe(23);
		expect(m3.mat[2][1]).toBe(-37);
		expect(m3.mat[2][2]).toBe(3);
	});

	it("mul prototype works with number as parameter", () => {
		let m3 = m1.mul(2);

		expect(m3.mat[0][0]).toBe(8);
		expect(m3.mat[0][1]).toBe(-4);
		expect(m3.mat[0][2]).toBe(6);
		expect(m3.mat[1][0]).toBe(2);
		expect(m3.mat[1][1]).toBe(4);
		expect(m3.mat[1][2]).toBe(-6);
		expect(m3.mat[2][0]).toBe(0);
		expect(m3.mat[2][1]).toBe(-2);
		expect(m3.mat[2][2]).toBe(14);
	});

	it("mul prototype throws TypeError when parameter is not a number or Mat3", () => {
		expect(() => new Mat3().mul("a")).toThrowError(TypeError);
		expect(() => new Mat3().mul()).toThrowError(TypeError);
	});

	it("transpose prototype works as expected", () => {
		let m3 = m1.transpose();

		expect(m3.mat[0][0]).toBe(4);
		expect(m3.mat[0][1]).toBe(1);
		expect(m3.mat[0][2]).toBe(0);
		expect(m3.mat[1][0]).toBe(-2);
		expect(m3.mat[1][1]).toBe(2);
		expect(m3.mat[1][2]).toBe(-1);
		expect(m3.mat[2][0]).toBe(3);
		expect(m3.mat[2][1]).toBe(-3);
		expect(m3.mat[2][2]).toBe(7);
	});

	it("det prototype works as expected", () => {
		let det = m1.det();
		expect(det).toBe(55);
	});

	it("inverse prototype works as expected", () => {
		let m3 = m1.inverse();

		expect(m3.mat[0][0]).toBeCloseTo(0.2, 5);
		expect(m3.mat[0][1]).toBeCloseTo(0.2, 5);
		expect(m3.mat[0][2]).toBeCloseTo(0, 5);
		expect(m3.mat[1][0]).toBeCloseTo(-0.12727272727272726, 5);
		expect(m3.mat[1][1]).toBeCloseTo(0.509090909090909, 5);
		expect(m3.mat[1][2]).toBeCloseTo(0.2727272727272727, 5);
		expect(m3.mat[2][0]).toBeCloseTo(-0.01818181818181818, 5);
		expect(m3.mat[2][1]).toBeCloseTo(0.07272727272727272, 5);
		expect(m3.mat[2][2]).toBeCloseTo(0.18181818181818182, 5);
	});

	it("c function returns itself", () => {
		expect(new Mat3().c() instanceof Mat3).toBe(true);
	});
});


describe("Mat3Identity", () => {
	let m = new Mat3Identity();

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(1);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(1);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat3).toBe(true);
		expect(typeof Mat3Identity.prototype.constructor).toBe("function");
		expect(typeof Mat3Identity.prototype.parent).toBe("function");
	});
});


describe("Mat3RotX", () => {
	let m = new Mat3RotX(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[1][2]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[2][2]).toBeCloseTo(-0.6669380616522619, 5);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat3RotX()).toThrowError(TypeError);
		expect(() => new Mat3RotX("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat3).toBe(true);
		expect(typeof Mat3RotX.prototype.constructor).toBe("function");
		expect(typeof Mat3RotX.prototype.parent).toBe("function");
	});
});


describe("Mat3RotY", () => {
	let m = new Mat3RotY(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(1);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[2][0]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBeCloseTo(-0.6669380616522619, 5);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat3RotY()).toThrowError(TypeError);
		expect(() => new Mat3RotY("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat3).toBe(true);
		expect(typeof Mat3RotY.prototype.constructor).toBe("function");
		expect(typeof Mat3RotY.prototype.parent).toBe("function");
	});
});


describe("Mat3RotZ", () => {
	let m = new Mat3RotZ(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[0][1]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[1][0]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[1][1]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(1);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat3RotZ()).toThrowError(TypeError);
		expect(() => new Mat3RotZ("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat3).toBe(true);
		expect(typeof Mat3RotZ.prototype.constructor).toBe("function");
		expect(typeof Mat3RotZ.prototype.parent).toBe("function");
	});
});


describe("Mat4", () => {
	let p1 = new Point3D(4, -2, 3, 4);
	let p2 = new Point3D(1, 2, -3, -1);
	let p3 = new Point3D(6, 3, 2, -3);
	let p4 = new Point3D(0, -1, 7, 1);
	let m1 = new Mat4(p1, p2, p3, p4);
	let p5 = new Point3D(1, -5, 1, 0);
	let p6 = new Point3D(-2, 2, 4, 1);
	let p7 = new Point3D(3, -5, 1, -2);
	let p8 = new Point3D(2, -4, 5, 5);
	let m2 = new Mat4(p5, p6, p7, p8);

	it("constructor works correctly with Point3D as parameters", () => {
		expect(m1.mat[0][0]).toBe(4);
		expect(m1.mat[0][1]).toBe(-2);
		expect(m1.mat[0][2]).toBe(3);
		expect(m1.mat[0][3]).toBe(4);
		expect(m1.mat[1][0]).toBe(1);
		expect(m1.mat[1][1]).toBe(2);
		expect(m1.mat[1][2]).toBe(-3);
		expect(m1.mat[1][3]).toBe(-1);
		expect(m1.mat[2][0]).toBe(6);
		expect(m1.mat[2][1]).toBe(3);
		expect(m1.mat[2][2]).toBe(2);
		expect(m1.mat[2][3]).toBe(-3);
		expect(m1.mat[3][0]).toBe(0);
		expect(m1.mat[3][1]).toBe(-1);
		expect(m1.mat[3][2]).toBe(7);
		expect(m1.mat[3][3]).toBe(1);
	});

	it("constructor works correctly with Mat4 as parameter", () => {
		let m2 = new Mat4(m1);

		expect(m2.mat[0][0]).toBe(4);
		expect(m2.mat[0][1]).toBe(-2);
		expect(m2.mat[0][2]).toBe(3);
		expect(m2.mat[0][3]).toBe(4);
		expect(m2.mat[1][0]).toBe(1);
		expect(m2.mat[1][1]).toBe(2);
		expect(m2.mat[1][2]).toBe(-3);
		expect(m2.mat[1][3]).toBe(-1);
		expect(m2.mat[2][0]).toBe(6);
		expect(m2.mat[2][1]).toBe(3);
		expect(m2.mat[2][2]).toBe(2);
		expect(m2.mat[2][3]).toBe(-3);
		expect(m2.mat[3][0]).toBe(0);
		expect(m2.mat[3][1]).toBe(-1);
		expect(m2.mat[3][2]).toBe(7);
		expect(m2.mat[3][3]).toBe(1);
	});

	it("constructor works correctly w/o parameters", () => {
		let m = new Mat4();

		expect(m.mat[0][0]).toBe(0);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(0);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(0);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(0);
	});

	it("add prototype works with Mat4 as parameter", () => {
		let m3 = m1.add(m2);

		expect(m3.mat[0][0]).toBe(5);
		expect(m3.mat[0][1]).toBe(-7);
		expect(m3.mat[0][2]).toBe(4);
		expect(m3.mat[0][3]).toBe(4);
		expect(m3.mat[1][0]).toBe(-1);
		expect(m3.mat[1][1]).toBe(4);
		expect(m3.mat[1][2]).toBe(1);
		expect(m3.mat[1][3]).toBe(0);
		expect(m3.mat[2][0]).toBe(9);
		expect(m3.mat[2][1]).toBe(-2);
		expect(m3.mat[2][2]).toBe(3);
		expect(m3.mat[2][3]).toBe(-5);
		expect(m3.mat[3][0]).toBe(2);
		expect(m3.mat[3][1]).toBe(-5);
		expect(m3.mat[3][2]).toBe(12);
		expect(m3.mat[3][3]).toBe(6);
	});

	it("add prototype throws TypeError when parameter is not a Mat4", () => {
		expect(() => new Mat4().add("a")).toThrowError(TypeError);
		expect(() => new Mat4().add()).toThrowError(TypeError);
	});

	it("mul prototype works with Mat4 as parameter", () => {
		let m3 = m1.mul(m2);

		expect(m3.mat[0][0]).toBe(25);
		expect(m3.mat[0][1]).toBe(-55);
		expect(m3.mat[0][2]).toBe(19);
		expect(m3.mat[0][3]).toBe(12);
		expect(m3.mat[1][0]).toBe(-14);
		expect(m3.mat[1][1]).toBe(18);
		expect(m3.mat[1][2]).toBe(1);
		expect(m3.mat[1][3]).toBe(3);
		expect(m3.mat[2][0]).toBe(0);
		expect(m3.mat[2][1]).toBe(-22);
		expect(m3.mat[2][2]).toBe(5);
		expect(m3.mat[2][3]).toBe(-16);
		expect(m3.mat[3][0]).toBe(25);
		expect(m3.mat[3][1]).toBe(-41);
		expect(m3.mat[3][2]).toBe(8);
		expect(m3.mat[3][3]).toBe(-10);
	});

	it("mul prototype works with number as parameter", () => {
		let m3 = m1.mul(2);

		expect(m3.mat[0][0]).toBe(8);
		expect(m3.mat[0][1]).toBe(-4);
		expect(m3.mat[0][2]).toBe(6);
		expect(m3.mat[0][3]).toBe(8);
		expect(m3.mat[1][0]).toBe(2);
		expect(m3.mat[1][1]).toBe(4);
		expect(m3.mat[1][2]).toBe(-6);
		expect(m3.mat[1][3]).toBe(-2);
		expect(m3.mat[2][0]).toBe(12);
		expect(m3.mat[2][1]).toBe(6);
		expect(m3.mat[2][2]).toBe(4);
		expect(m3.mat[2][3]).toBe(-6);
		expect(m3.mat[3][0]).toBe(0);
		expect(m3.mat[3][1]).toBe(-2);
		expect(m3.mat[3][2]).toBe(14);
		expect(m3.mat[3][3]).toBe(2);
	});

	it("mul prototype throws TypeError when parameter is not a number or Mat4", () => {
		expect(() => new Mat4().mul("a")).toThrowError(TypeError);
		expect(() => new Mat4().mul()).toThrowError(TypeError);
	});

	it("transpose prototype works as expected", () => {
		let m3 = m1.transpose();

		expect(m3.mat[0][0]).toBe(4);
		expect(m3.mat[0][1]).toBe(1);
		expect(m3.mat[0][2]).toBe(6);
		expect(m3.mat[0][3]).toBe(0);
		expect(m3.mat[1][0]).toBe(-2);
		expect(m3.mat[1][1]).toBe(2);
		expect(m3.mat[1][2]).toBe(3);
		expect(m3.mat[1][3]).toBe(-1);
		expect(m3.mat[2][0]).toBe(3);
		expect(m3.mat[2][1]).toBe(-3);
		expect(m3.mat[2][2]).toBe(2);
		expect(m3.mat[2][3]).toBe(7);
		expect(m3.mat[3][0]).toBe(4);
		expect(m3.mat[3][1]).toBe(-1);
		expect(m3.mat[3][2]).toBe(-3);
		expect(m3.mat[3][3]).toBe(1);
	});

	it("toMat3 prototype works as expected", () => {
		let m2 = m1.toMat3();
		expect(m2.mat[0][0]).toBe(4);
		expect(m2.mat[0][1]).toBe(-2);
		expect(m2.mat[0][2]).toBe(3);
		expect(m2.mat[1][0]).toBe(1);
		expect(m2.mat[1][1]).toBe(2);
		expect(m2.mat[1][2]).toBe(-3);
		expect(m2.mat[2][0]).toBe(6);
		expect(m2.mat[2][1]).toBe(3);
		expect(m2.mat[2][2]).toBe(2);
	});

	it("c function returns itself", () => {
		expect(new Mat3().c() instanceof Mat3).toBe(true);
	});
});


describe("Mat4Identity", () => {
	let m = new Mat4Identity();

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(1);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(1);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4Identity.prototype.constructor).toBe("function");
		expect(typeof Mat4Identity.prototype.parent).toBe("function");
	});
});


describe("Mat4OrthoRH", () => {
	let m = new Mat4OrthoRH(800, 600, 0.1, 100);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(0.0025);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBeCloseTo(0.0033333333333333335, 5);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBeCloseTo(-0.01001001001001001, 5);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(-1);
		expect(m.mat[3][1]).toBe(-1);
		expect(m.mat[3][2]).toBeCloseTo(-0.001001001001001001, 5);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4OrthoRH(1, 2, 3)).toThrowError(TypeError);
		expect(() => new Mat4OrthoRH(1, 2, 3, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4OrthoRH.prototype.constructor).toBe("function");
		expect(typeof Mat4OrthoRH.prototype.parent).toBe("function");
	});
});


describe("Mat4PerspRH", () => {
	let m = new Mat4PerspRH(60, 800/600, 0.1, 100);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.20815993621554563, 5);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBeCloseTo(-0.15611995216165922, 5);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBeCloseTo(-1.0010010010010009, 5);
		expect(m.mat[2][3]).toBe(-1);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBeCloseTo(-0.10010010010010009, 5);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4PerspRH(1, 2, 3)).toThrowError(TypeError);
		expect(() => new Mat4PerspRH(1, 2, 3, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4PerspRH.prototype.constructor).toBe("function");
		expect(typeof Mat4PerspRH.prototype.parent).toBe("function");
	});
});


describe("Mat4RotX", () => {
	let m = new Mat4RotX(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[1][2]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[2][2]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat4RotX()).toThrowError(TypeError);
		expect(() => new Mat4RotX("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4RotX.prototype.constructor).toBe("function");
		expect(typeof Mat4RotX.prototype.parent).toBe("function");
	});
});


describe("Mat4RotY", () => {
	let m = new Mat4RotY(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[0][2]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(1);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat4RotY()).toThrowError(TypeError);
		expect(() => new Mat4RotY("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4RotY.prototype.constructor).toBe("function");
		expect(typeof Mat4RotY.prototype.parent).toBe("function");
	});
});


describe("Mat4RotZ", () => {
	let m = new Mat4RotZ(40);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[0][1]).toBeCloseTo(0.7451131604793488, 5);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBeCloseTo(-0.7451131604793488, 5);
		expect(m.mat[1][1]).toBeCloseTo(-0.6669380616522619, 5);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(1);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameter is not a number", () => {
		expect(() => new Mat4RotZ()).toThrowError(TypeError);
		expect(() => new Mat4RotZ("a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4RotZ.prototype.constructor).toBe("function");
		expect(typeof Mat4RotZ.prototype.parent).toBe("function");
	});
});


describe("Mat4RotXYZ", () => {
	let m = new Mat4RotXYZ(40, 50, 110);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(-0.9640211466051968, 5);
		expect(m.mat[0][1]).toBeCloseTo(-0.042692681361605984, 5);
		expect(m.mat[0][2]).toBeCloseTo(0.26237485370392877, 5);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBeCloseTo(0.16580040053409972, 5);
		expect(m.mat[1][1]).toBeCloseTo(0.6749344021795681, 5);
		expect(m.mat[1][2]).toBeCloseTo(0.7190088872449639, 5);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBeCloseTo(-0.20778223235092494, 5);
		expect(m.mat[2][1]).toBeCloseTo(0.736641627735404, 5);
		expect(m.mat[2][2]).toBeCloseTo(-0.6435725726028114, 5);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4RotXYZ(1, 2)).toThrowError(TypeError);
		expect(() => new Mat4RotXYZ(1, 2, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4RotXYZ.prototype.constructor).toBe("function");
		expect(typeof Mat4RotXYZ.prototype.parent).toBe("function");
	});
});


describe("Mat4Scale", () => {
	let m = new Mat4Scale(1.5, 2, -1);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1.5);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(2);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(-1);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(0);
		expect(m.mat[3][1]).toBe(0);
		expect(m.mat[3][2]).toBe(0);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4Scale(1, 2)).toThrowError(TypeError);
		expect(() => new Mat4Scale(1, 2, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4Scale.prototype.constructor).toBe("function");
		expect(typeof Mat4Scale.prototype.parent).toBe("function");
	});
});


describe("Mat4Transl", () => {
	let m = new Mat4Transl(1.5, 2, -1);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBe(1);
		expect(m.mat[0][1]).toBe(0);
		expect(m.mat[0][2]).toBe(0);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBe(0);
		expect(m.mat[1][1]).toBe(1);
		expect(m.mat[1][2]).toBe(0);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBe(0);
		expect(m.mat[2][1]).toBe(0);
		expect(m.mat[2][2]).toBe(1);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBe(1.5);
		expect(m.mat[3][1]).toBe(2);
		expect(m.mat[3][2]).toBe(-1);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4Transl(1, 2)).toThrowError(TypeError);
		expect(() => new Mat4Transl(1, 2, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4Transl.prototype.constructor).toBe("function");
		expect(typeof Mat4Transl.prototype.parent).toBe("function");
	});
});


describe("Mat4ViewRH", () => {
	let v1 = new Vec3D(1, 2, 3);
	let v2 = new Vec3D(-4, 2, -3);
	let v3 = new Vec3D(1.5, 7, -1);
	let m = new Mat4ViewRH(v1, v2, v3);

	it("constructor works", () => {
		expect(m.mat[0][0]).toBeCloseTo(0.5088423789474245, 5);
		expect(m.mat[0][1]).toBeCloseTo(0.4351497391198703, 5);
		expect(m.mat[0][2]).toBeCloseTo(0.7427813527082074, 5);
		expect(m.mat[0][3]).toBe(0);
		expect(m.mat[1][0]).toBeCloseTo(-0.22764001163437408, 5);
		expect(m.mat[1][1]).toBeCloseTo(0.9001383174936746, 5);
		expect(m.mat[1][2]).toBeCloseTo(-0.3713906763541037, 5);
		expect(m.mat[1][3]).toBe(0);
		expect(m.mat[2][0]).toBeCloseTo(-0.8302165130194821, 5);
		expect(m.mat[2][1]).toBeCloseTo(0.01989255950262267, 5);
		expect(m.mat[2][2]).toBeCloseTo(0.5570860145311556, 5);
		expect(m.mat[2][3]).toBe(0);
		expect(m.mat[3][0]).toBeCloseTo(2.43708718337977, 5);
		expect(m.mat[3][1]).toBeCloseTo(-2.295104052615087 , 5);
		expect(m.mat[3][2]).toBeCloseTo(-1.6712580435934667, 5);
		expect(m.mat[3][3]).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Mat4ViewRH(v1, v2)).toThrowError(TypeError);
		expect(() => new Mat4ViewRH(v1, v2, "a")).toThrowError(TypeError);
	});

	it("prototype is set correctly", () => {
		expect(m instanceof Mat4).toBe(true);
		expect(typeof Mat4ViewRH.prototype.constructor).toBe("function");
		expect(typeof Mat4ViewRH.prototype.parent).toBe("function");
	});
});


describe("Camera", () => {

	it("constructor sets all variables correctly", () => {
		let cam = new Camera();
		expect(cam.azimuth).toBe(0);
		expect(cam.zenith).toBe(0);
		expect(cam.radius).toBe(1);
		expect(cam.xy).toBe(true);
		expect(cam.pos instanceof Vec3D).toBe(true);
		expect(cam.firstPerson).toBe(true);
		expect(cam.eye instanceof Vec3D).toBe(true);
		expect(cam.eyeVector instanceof Vec3D).toBe(true);
		expect(cam.view instanceof Mat4).toBe(true);
	});

	it("computeMatrix prototype works correctly with firstPerson==true", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.setZenith(-0.1);
		cam.setRadius(1.1);
		cam.setPosition(new Vec3D(1, 5, -4.5));
		cam.setFirstPerson(true);

		expect(cam.eyeVector.x).toBeCloseTo(-0.4770304078518429, 5);
		expect(cam.eyeVector.y).toBeCloseTo(0.8731983044562817, 5);
		expect(cam.eyeVector.z).toBeCloseTo(-0.09983341664682815, 5);

		expect(cam.eye.x).toBe(1);
		expect(cam.eye.y).toBe(5);
		expect(cam.eye.z).toBe(-4.5);

		expect(cam.view.mat[0][0]).toBeCloseTo(0.8775825618903728, 5);
		expect(cam.view.mat[0][1]).toBeCloseTo(-0.04786268954660339, 5);
		expect(cam.view.mat[0][2]).toBeCloseTo(0.47703040785184286, 5);
		expect(cam.view.mat[0][3]).toBe(0);

		expect(cam.view.mat[1][0]).toBeCloseTo(0.47942553860420295, 5);
		expect(cam.view.mat[1][1]).toBeCloseTo(0.08761206554319244, 5);
		expect(cam.view.mat[1][2]).toBeCloseTo(-0.8731983044562817, 5);
		expect(cam.view.mat[1][3]).toBe(0);

		expect(cam.view.mat[2][0]).toBeCloseTo(0, 5);
		expect(cam.view.mat[2][1]).toBeCloseTo(0.9950041652780257, 5);
		expect(cam.view.mat[2][2]).toBeCloseTo(0.09983341664682815, 5);
		expect(cam.view.mat[2][3]).toBe(0);

		expect(cam.view.mat[3][0]).toBeCloseTo(-3.2747102549113873, 5);
		expect(cam.view.mat[3][1]).toBeCloseTo(4.087321105581757, 5);
		expect(cam.view.mat[3][2]).toBeCloseTo(4.338211489340292, 5);
		expect(cam.view.mat[3][3]).toBe(1);
	});

	it("computeMatrix prototype works correctly with firstPerson==false", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.setZenith(-0.1);
		cam.setRadius(1.1);
		cam.setPosition(new Vec3D(1, 5, -4.5));
		cam.setFirstPerson(false);

		expect(cam.eyeVector.x).toBeCloseTo(-0.4770304078518429, 5);
		expect(cam.eyeVector.y).toBeCloseTo(0.8731983044562817, 5);
		expect(cam.eyeVector.z).toBeCloseTo(-0.09983341664682815, 5);

		expect(cam.eye.x).toBeCloseTo(1.5247334486370272, 5);
		expect(cam.eye.y).toBeCloseTo(4.03948186509809, 5);
		expect(cam.eye.z).toBeCloseTo(-4.390183241688489, 5);

		expect(cam.view.mat[0][0]).toBeCloseTo(0.8775825618903728, 5);
		expect(cam.view.mat[0][1]).toBeCloseTo(-0.04786268954660339, 5);
		expect(cam.view.mat[0][2]).toBeCloseTo(0.47703040785184286, 5);
		expect(cam.view.mat[0][3]).toBe(0);

		expect(cam.view.mat[1][0]).toBeCloseTo(0.47942553860420295, 5);
		expect(cam.view.mat[1][1]).toBeCloseTo(0.08761206554319244, 5);
		expect(cam.view.mat[1][2]).toBeCloseTo(-0.8731983044562817, 5);
		expect(cam.view.mat[1][3]).toBe(0);

		expect(cam.view.mat[2][0]).toBeCloseTo(0, 5);
		expect(cam.view.mat[2][1]).toBeCloseTo(0.9950041652780257, 5);
		expect(cam.view.mat[2][2]).toBeCloseTo(0.09983341664682815, 5);
		expect(cam.view.mat[2][3]).toBe(0);

		expect(cam.view.mat[3][0]).toBeCloseTo(-3.2747102549113873, 5);
		expect(cam.view.mat[3][1]).toBeCloseTo(4.087321105581757, 5);
		expect(cam.view.mat[3][2]).toBeCloseTo(3.238211489340292, 5);
		expect(cam.view.mat[3][3]).toBe(1);
	});

	it("addAzimuth prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.addAzimuth(2);
		cam.addAzimuth(5);
		expect(cam.azimuth).toBe(7);

		expect(cam.eyeVector.x).toBeCloseTo(-0.6569865987187891, 5);
	});

	it("addAzimuth prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().addAzimuth("a")).toThrowError(TypeError);
		expect(() => new Camera().addAzimuth()).toThrowError(TypeError);
	});

	it("setAzimuth prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setAzimuth(2);
		expect(cam.azimuth).toBe(2);
		cam.setAzimuth(5);
		expect(cam.azimuth).toBe(5);

		expect(cam.eyeVector.x).toBeCloseTo(0.9589242746631385 , 5);
	});

	it("setAzimuth prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().setAzimuth("a")).toThrowError(TypeError);
		expect(() => new Camera().setAzimuth()).toThrowError(TypeError);
	});

	it("addZenith prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.addZenith(2);
		expect(cam.zenith).toBe(0);
		expect(cam.eyeVector.y).toBe(1);
		cam.addZenith(Math.PI / 2);
		expect(cam.zenith).toBe(Math.PI / 2);
	});

	it("addZenith prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().addZenith("a")).toThrowError(TypeError);
		expect(() => new Camera().addZenith()).toThrowError(TypeError);
	});

	it("setZenith prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setZenith(2);
		expect(cam.zenith).toBe(2);
		cam.setZenith(5);
		expect(cam.zenith).toBe(5);

		expect(cam.eyeVector.y).toBeCloseTo(0.28366218546322625, 5);
	});

	it("setZenith prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().setZenith("a")).toThrowError(TypeError);
		expect(() => new Camera().setZenith()).toThrowError(TypeError);
	});

	it("addRadius prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.addRadius(2);
		cam.addRadius(5);
		expect(cam.radius).toBe(8);
		cam.addRadius(-8);
		expect(cam.radius).toBe(0.1);

		cam.setFirstPerson(false);
		cam.addRadius(0.3);
		expect(cam.eye.y).toBe(-0.4);
	});

	it("addRadius prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().addRadius("a")).toThrowError(TypeError);
		expect(() => new Camera().addRadius()).toThrowError(TypeError);
	});

	it("setRadius prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setRadius(2);
		expect(cam.radius).toBe(2);
		cam.setRadius(5);
		expect(cam.radius).toBe(5);

		cam.setFirstPerson(false);
		cam.setRadius(7);
		expect(cam.eye.y).toBe(-7);
	});

	it("setRadius prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().setRadius("a")).toThrowError(TypeError);
		expect(() => new Camera().setRadius()).toThrowError(TypeError);
	});

	it("mulRadius prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.mulRadius(2);
		cam.mulRadius(5);
		expect(cam.radius).toBe(10);
		cam.mulRadius(0);
		expect(cam.radius).toBe(0.1);

		cam.setFirstPerson(false);
		cam.mulRadius(1.1);
		expect(cam.eye.y).toBeCloseTo(-0.11, 5);
	});

	it("mulRadius prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().mulRadius("a")).toThrowError(TypeError);
		expect(() => new Camera().mulRadius()).toThrowError(TypeError);
	});

	it("forward prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.forward(2);
		expect(cam.pos.x).toBeCloseTo(-0.9588510772084061, 5);
		expect(cam.pos.y).toBeCloseTo(1.7551651237807453, 5);
		expect(cam.pos.z).toBe(0);

		cam.setFirstPerson(false);
		cam.forward(1.5);
		expect(cam.pos.x).toBeCloseTo(-1.6779893851147107, 5);
		expect(cam.pos.y).toBeCloseTo(3.071538966616304, 5);
		expect(cam.pos.z).toBe(0);

		cam.setFirstPerson(true);
		cam.forward(1.5);
		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("forward prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().forward("a")).toThrowError(TypeError);
		expect(() => new Camera().forward()).toThrowError(TypeError);
	});

	it("backward prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.backward(2);
		expect(cam.pos.x).toBeCloseTo(0.9588510772084061, 5);
		expect(cam.pos.y).toBeCloseTo(-1.7551651237807453, 5);
		expect(cam.pos.z).toBe(0);

		cam.setFirstPerson(false);
		cam.backward(1.5);
		expect(cam.pos.x).toBeCloseTo(1.6779893851147107, 5);
		expect(cam.pos.y).toBeCloseTo(-3.071538966616304, 5);
		expect(cam.pos.z).toBe(0);

		cam.setFirstPerson(true);
		cam.backward(1.5);
		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("backward prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().backward("a")).toThrowError(TypeError);
		expect(() => new Camera().backward()).toThrowError(TypeError);
	});

	it("right prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.right(2);
		expect(cam.pos.x).toBeCloseTo(1.7551651237807455, 5);
		expect(cam.pos.y).toBeCloseTo(0.958851077208406, 5);
		expect(cam.pos.z).toBe(0);

		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("right prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().right("a")).toThrowError(TypeError);
		expect(() => new Camera().right()).toThrowError(TypeError);
	});

	it("left prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.setAzimuth(0.5);
		cam.left(2);
		expect(cam.pos.x).toBeCloseTo(-1.7551651237807455, 5);
		expect(cam.pos.y).toBeCloseTo(-0.958851077208406, 5);
		expect(cam.pos.z).toBe(0);

		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("left prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().left("a")).toThrowError(TypeError);
		expect(() => new Camera().left()).toThrowError(TypeError);
	});

	it("down prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.down(2);
		expect(cam.pos.x).toBe(0);
		expect(cam.pos.y).toBe(0);
		expect(cam.pos.z).toBe(-2);

		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("down prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().down("a")).toThrowError(TypeError);
		expect(() => new Camera().down()).toThrowError(TypeError);
	});

	it("up prototype works with number as a parameter", () => {
		let cam = new Camera();
		cam.up(2);
		expect(cam.pos.x).toBe(0);
		expect(cam.pos.y).toBe(0);
		expect(cam.pos.z).toBe(2);

		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("up prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Camera().up("a")).toThrowError(TypeError);
		expect(() => new Camera().up()).toThrowError(TypeError);
	});

	it("move prototype works with Vec3D as a parameter", () => {
		let cam = new Camera();
		cam.move(new Vec3D(1, 2, 3));
		cam.move(new Vec3D(1.5, -2, -1));
		expect(cam.pos.x).toBe(2.5);
		expect(cam.pos.y).toBe(0);
		expect(cam.pos.z).toBe(2);

		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("move prototype throws TypeError when parameter is not a Vec3D", () => {
		expect(() => new Camera().move("a")).toThrowError(TypeError);
		expect(() => new Camera().move()).toThrowError(TypeError);
	});

	it("setPosition prototype works with Vec3D as a parameter", () => {
		let cam = new Camera();
		cam.setPosition(new Vec3D(1.5, -2, -1));
		expect(cam.pos.x).toBe(1.5);
		expect(cam.pos.y).toBe(-2);
		expect(cam.pos.z).toBe(-1);

		expect(cam.pos.x).toBe(cam.eye.x);
		expect(cam.pos.y).toBe(cam.eye.y);
		expect(cam.pos.z).toBe(cam.eye.z);
	});

	it("setPosition prototype throws TypeError when parameter is not a Vec3D", () => {
		expect(() => new Camera().setPosition("a")).toThrowError(TypeError);
		expect(() => new Camera().setPosition()).toThrowError(TypeError);
	});

	it("setFirstPerson prototype works with boolean as a parameter", () => {
		let cam = new Camera();
		cam.setFirstPerson(false);
		expect(cam.firstPerson).toBe(false);
		cam.setFirstPerson(true);
		expect(cam.firstPerson).toBe(true);
	});

	it("setFirstPerson prototype throws TypeError when parameter is not a boolean", () => {
		expect(() => new Camera().setFirstPerson("a")).toThrowError(TypeError);
		expect(() => new Camera().setFirstPerson()).toThrowError(TypeError);
	});
});


describe("Col", () => {

	it("constructor can make a copy of another Col", () => {
		let v = new Col(new Col(0.1, 0.2, 0.3, 1));
		expect(v.r).toBe(0.1);
		expect(v.g).toBe(0.2);
		expect(v.b).toBe(0.3);
		expect(v.a).toBe(1);
	});

	it("constructor can make a copy of a Point3D", () => {
		let v = new Col(new Point3D(5, 7, 8, 9));
		expect(v.r).toBe(5);
		expect(v.g).toBe(7);
		expect(v.b).toBe(8);
		expect(v.a).toBe(9);
	});

	it("constructor with integer parameters works", () => {
		let v1 = new Col(50, 100, 150, 200);
		expect(v1.r).toBe(50/255);
		expect(v1.g).toBe(100/255);
		expect(v1.b).toBe(150/255);
		expect(v1.a).toBe(200/255);

		let v2 = new Col(50, 100, 150);
		expect(v2.r).toBe(50/255);
		expect(v2.g).toBe(100/255);
		expect(v2.b).toBe(150/255);
		expect(v2.a).toBe(1);
	});

	it("constructor with float parameters works", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 0.4);
		expect(v1.r).toBe(0.1);
		expect(v1.g).toBe(0.2);
		expect(v1.b).toBe(0.3);
		expect(v1.a).toBe(0.4);

		let v2 = new Col(0.1, 0.2, 0.3);
		expect(v2.r).toBe(0.1);
		expect(v2.g).toBe(0.2);
		expect(v2.b).toBe(0.3);
		expect(v2.a).toBe(1);
	});

	it("constructor throws TypeError when parameters are wrong", () => {
		expect(() => new Col("a")).toThrowError(TypeError);
		expect(() => new Col(20, 20)).toThrowError(TypeError);
		expect(() => new Col()).toThrowError(TypeError);
	});

	it("addna prototype works as expected", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 100);
		let v2 = new Col(0.2, 0.4, 0.6, 100);
		let v3 = v1.addna(v2);
		expect(v3.r).toBeCloseTo(0.3, 5);
		expect(v3.g).toBeCloseTo(0.6, 5);
		expect(v3.b).toBeCloseTo(0.9, 5);
		expect(v3.a).toBeCloseTo(1, 5);
	});

	it("addna prototype throws TypeError when parameter is not Col", () => {
		expect(() => new Col(1, 2, 3, 4).addna()).toThrowError(TypeError);
		expect(() => new Col(1, 2, 3, 4).addna(1)).toThrowError(TypeError);
	});

	it("mulna prototype works as expected", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 100);
		let v3 = v1.mulna(2);
		expect(v3.r).toBeCloseTo(0.2, 5);
		expect(v3.g).toBeCloseTo(0.4, 5);
		expect(v3.b).toBeCloseTo(0.6, 5);
		expect(v3.a).toBeCloseTo(1, 5);
	});

	it("mulna prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Col(1, 2, 3, 4).mulna()).toThrowError(TypeError);
		expect(() => new Col(1, 2, 3, 4).mulna("a")).toThrowError(TypeError);
	});

	it("add prototype works as expected", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 100);
		let v2 = new Col(0.2, 0.4, 0.6, 100);
		let v3 = v1.add(v2);
		expect(v3.r).toBeCloseTo(0.3, 5);
		expect(v3.g).toBeCloseTo(0.6, 5);
		expect(v3.b).toBeCloseTo(0.9, 5);
		expect(v3.a).toBeCloseTo(200, 5);
	});

	it("add prototype throws TypeError when parameter is not Col", () => {
		expect(() => new Col(1, 2, 3, 4).add()).toThrowError(TypeError);
		expect(() => new Col(1, 2, 3, 4).add(1)).toThrowError(TypeError);
	});

	it("mul prototype works as expected", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 100);
		let v2 = new Col(0.2, 0.4, 0.6, 100);
		let v3 = v1.mul(v2);
		expect(v3.r).toBeCloseTo(0.02, 5);
		expect(v3.g).toBeCloseTo(0.08, 5);
		expect(v3.b).toBeCloseTo(0.18, 5);
		expect(v3.a).toBeCloseTo(10000, 5);

		let v4 = v2.mul(2.5);
		expect(v4.r).toBeCloseTo(0.5, 5);
		expect(v4.g).toBeCloseTo(1.0, 5);
		expect(v4.b).toBeCloseTo(1.5, 5);
		expect(v4.a).toBeCloseTo(250, 5);
	});

	it("mul prototype throws TypeError when parameter is not Col or number", () => {
		expect(() => new Col(1, 2, 3, 4).mul()).toThrowError(TypeError);
		expect(() => new Col(1, 2, 3, 4).mul("a")).toThrowError(TypeError);
	});

	it("gamma prototype works as expected", () => {
		let v1 = new Col(0.1, 0.2, 0.3, 100);
		let v2 = v1.gamma(3);
		expect(v2.r).toBeCloseTo(0.001, 5);
		expect(v2.g).toBeCloseTo(0.008, 5);
		expect(v2.b).toBeCloseTo(0.027, 5);
		expect(v2.a).toBeCloseTo(100, 5);
	});

	it("gamma prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Col(1, 2, 3, 4).gamma()).toThrowError(TypeError);
		expect(() => new Col(1, 2, 3, 4).gamma("a")).toThrowError(TypeError);
	});

	it("saturate prototype works as expected", () => {
		let v1 = new Col(1000, -1000, 10, 100);
		let v2 = v1.saturate();
		expect(v2.r).toBe(1);
		expect(v2.g).toBe(0);
		expect(v2.b).toBe(10/255);
		expect(v2.a).toBe(100/255);
	});

	it("getRGB prototype works as expected", () => {
		let v1 = new Col(10, 20, 30, 40);
		let rgb1 = v1.getRGB();
		expect(rgb1).toBe(660510);

		let v2 = new Col(0.5, 0, 1, 0);
		let rgb2 = v2.getRGB();
		expect(rgb2).toBe(8323327);
	});

	it("getARGB prototype works as expected", () => {
		let v1 = new Col(10, 20, 30, 40);
		let rgb1 = v1.getARGB();
		expect(rgb1).toBe(671749150);

		let v2 = new Col(0.5, 0, 1, 0.1);
		let rgb2 = v2.getARGB();
		expect(rgb2).toBe(427753727);
	});

	it("c function returns itself", () => {
		expect(new Col(1, 1, 1, 1).c() instanceof Col).toBe(true);
	});
});


describe("Kubika", () => {
	let p1 = new Point3D(1, 2, 3, 4);
	let p2 = new Point3D(-4, -2, 0, 7);
	let p3 = new Point3D(5, -1, 4, -1);
	let p4 = new Point3D(0, 1.5, 4.5, -3);

	it("constructor works with all types of Kubika", () => {
		let k1 = new Kubika(0);
		expect(k1.bm.mat[0][0]).toBe(-1);
		expect(k1.bm.mat[2][2]).toBe(0);
		expect(k1.type).toBe(0);

		let k2 = new Kubika(1);
		expect(k2.bm.mat[0][0]).toBe(2);
		expect(k2.bm.mat[2][2]).toBe(1);
		expect(k2.type).toBe(1);

		let k3 = new Kubika(2);
		expect(k3.bm.mat[0][0]).toBeCloseTo(-1/6, 5);
		expect(k3.bm.mat[2][2]).toBeCloseTo(3/6, 5);
		expect(k3.type).toBe(2);

		let k4 = new Kubika();
		expect(k4.bm.mat[0][0]).toBe(-1);
		expect(k4.bm.mat[2][2]).toBe(0);
		expect(k4.type).toBe(0);
	});

	it("init prototype works as expected", () => {
		let k1 = new Kubika(0);
		k1.init(p1, p2, p3, p4);
		expect(k1.rb instanceof Mat4).toBe(true);
		expect(k1.rb.mat[0][0]).toBe(-28);
		expect(k1.rb.mat[2][1]).toBe(-12);
		expect(k1.rb.mat[1][3]).toBe(-33);

		let k2 = new Kubika(1);
		k2.init(p1, p2, p3, p4);
		expect(k2.rb instanceof Mat4).toBe(true);
		expect(k2.rb.mat[0][0]).toBe(-8);
		expect(k2.rb.mat[2][1]).toBe(-4);
		expect(k2.rb.mat[1][3]).toBe(-25);

		let k3 = new Kubika(2);
		k3.init(p1, p2, p3, p4);
		expect(k3.rb instanceof Mat4).toBe(true);
		expect(k3.rb.mat[0][0]).toBeCloseTo(-28/6, 5);
		expect(k3.rb.mat[2][1]).toBeCloseTo(-9/6, 5);
		expect(k3.rb.mat[1][3]).toBeCloseTo(-33/6, 5);
	});

	it("init prototype throws TypeError when parameters are not Point3D", () => {
		expect(() => new Kubika(0).init()).toThrowError(TypeError);
		expect(() => new Kubika(0).init("a", "a", "a", "a")).toThrowError(TypeError);
	});

	it("compute prototype works as expected", () => {
		let k1 = new Kubika(0);
		k1.init(p1, p2, p3, p4);
		let pp1 = k1.compute(-1);
		let pp2 = k1.compute(0);
		expect(pp1.x).toBe(pp2.x);
		expect(pp1.y).toBe(pp2.y);
		expect(pp1.z).toBe(pp2.z);
		expect(pp1.w).toBe(pp2.w);

		let pp3 = k1.compute(2);
		let pp4 = k1.compute(1);
		expect(pp3.x).toBe(pp4.x);
		expect(pp3.y).toBe(pp4.y);
		expect(pp3.z).toBe(pp4.z);
		expect(pp3.w).toBe(pp4.w);

		let pp5 = k1.compute(0.4);
		expect(pp5.x).toBeCloseTo(-0.072, 5);
		expect(pp5.y).toBeCloseTo(-0.624, 5);
		expect(pp5.z).toBeCloseTo(2.088, 5);
		expect(pp5.w).toBe(1);
	});

	it("compute prototype throws TypeError when parameter is not a number", () => {
		expect(() => new Kubika(0).compute()).toThrowError(TypeError);
		expect(() => new Kubika(0).compute("a")).toThrowError(TypeError);
	});
});


describe("Bikubika", () => {
	let p11 = new Point3D(1, 2, -3, 4);
	let p12 = new Point3D(-4, -2, 0, 7);
	let p13 = new Point3D(5, -1, 4, -1);
	let p14 = new Point3D(2, 5.5, 4.5, -3);
	let p21 = new Point3D(1, 2, 3, 4);
	let p22 = new Point3D(-4, -2, 1, 7);
	let p23 = new Point3D(5, -1, 4, -1);
	let p24 = new Point3D(0, 1.5, 4.5, -3);
	let p31 = new Point3D(1, 3, 3, 4);
	let p32 = new Point3D(0, -2, 0, 7);
	let p33 = new Point3D(5, -1, 4, -1);
	let p34 = new Point3D(3, 1.5, 4.5, -3);
	let p41 = new Point3D(1, 4, 3, 7);
	let p42 = new Point3D(-4, -2, 0, 7);
	let p43 = new Point3D(5, -5, 4, -1);
	let p44 = new Point3D(-2, 1.5, 0.5, -3);

	it("constructor sets all variables correctly", () => {
		let b1 = new Bikubika(0);
		expect(b1.u1).not.toBe("undefined");
		expect(b1.u2).not.toBe("undefined");
		expect(b1.u3).not.toBe("undefined");
		expect(b1.u4).not.toBe("undefined");
		expect(b1.k1 instanceof Kubika).toBe(true);
		expect(b1.k2 instanceof Kubika).toBe(true);
		expect(b1.k3 instanceof Kubika).toBe(true);
		expect(b1.k4 instanceof Kubika).toBe(true);
		expect(b1.k5 instanceof Kubika).toBe(true);
	});

	it("init prototype works as expected", () => {
		let b1 = new Bikubika(2);
		b1.init(p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43, p44);
		expect(b1.k1.rb.mat[0][0]).toBeCloseTo(-26/6, 5);
	});

	it("compute prototype works as expected", () => {
		let b1 = new Bikubika(2);
		b1.init(p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43, p44);

		let pp1 = b1.compute(-1, 0.5);
		let pp2 = b1.compute(0, 0.5);
		expect(pp1.x).toBe(pp2.x);
		expect(pp1.y).toBe(pp2.y);
		expect(pp1.z).toBe(pp2.z);
		expect(pp1.w).toBe(pp2.w);

		let pp3 = b1.compute(0.5, 2);
		let pp4 = b1.compute(0.5, 1);
		expect(pp3.x).toBe(pp4.x);
		expect(pp3.y).toBe(pp4.y);
		expect(pp3.z).toBe(pp4.z);
		expect(pp3.w).toBe(pp4.w);

		let pp5 = b1.compute(0.4, 0.7);
		expect(pp5.x).toBeCloseTo(1.2440408888888888, 5);
		expect(pp5.y).toBeCloseTo(-1.4732664444444443, 5);
		expect(pp5.z).toBeCloseTo(1.9988013333333332, 5);
		expect(pp5.w).toBe(1);
	});

	it("compute prototype throws TypeError when parameters are not numbers", () => {
		expect(() => new Kubika(0).compute()).toThrowError(TypeError);
		expect(() => new Kubika(0).compute("a", 1)).toThrowError(TypeError);
	});
});
