"use strict";

/**
 * @description This file contains set of objects and functions that allow to use vector and matrix operations in JavaScript.
 *  Then it also contains objects for work with cubics, bicubics, quaternions and camera.
 *
 *  It is based on https://gitlab.com/honza.vanek/transforms
 *  These functions are used in course of computer graphics (PGRF 1-3) at FIM UHK.
 *  (Faculty of Informatics and Management at University of Hradec Kralove)
 *
 * @version 1.2.1 (2020/03/02)
 * Backwards compatible with 1.1
 */

/**
 * Create square matrix with given size containing zeros.
 * @param {number} x matrix size (e.g. 3 for 3×3 matrix), if it is not filled then it is assumed 4×4
 * @throws {TypeError} If x is set and it is not a number
 * @constructor
 */
const ZeroArray = function(x) {
	const length = (typeof x !== "undefined") ? x : 4;
	if (typeof length === "number") {
		const mat = [];
		for (let i = 0; i < length; i++) {
			mat[i] = [];
			for (let j = 0; j < length; j++) {
				mat[i][j] = 0.0;
			}
		}
		return mat;
	} else {
		throw new TypeError("ZeroArray: Invalid parameter: must be a number");
	}
};

/**
 * Object for working with 1D vectors
 * @param {Vec1D,number} x number or Vec1D, if it is not filled then it is assumed 0
 * @constructor
 */
const Vec1D = function(x) {
	if (x instanceof Vec1D) {
		this.x = x.x;
	} else {
		if (typeof x !== "undefined") {
			if (typeof x === "number") {
				this.x = x;
			} else {
				throw new TypeError("Vec1D, constructor: Invalid parameter: must be a Vec1D or a number");
			}
		} else {
			this.x = 0.0;
		}
	}
};

/**
 * Add another vector to this vector
 * @param {Vec1D} v vector to be added
 * @return {Vec1D}  new Vec1D instance
 * @throws {TypeError} If v is not a Vec1D
 */
Vec1D.prototype.add = function(v) {
	if (v instanceof Vec1D) {
		return new Vec1D(this.x + v.x);
	} else {
		throw new TypeError("Vec1D, add: Invalid parameter: must be a Vec1D");
	}
};

/**
 * Multiplication by scalar
 * @param  {number} m scalar
 * @return {Vec1D}    new Vec1D instance
 * @throws {TypeError} If m is not a number
 */
Vec1D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec1D(this.x * m);
	} else {
		throw new TypeError("Vec1D, mul: Invalid parameter: must be a number");
	}
};

/**
 * Print value to the console
 * @return {Vec1D} reference to this instance for calls chaining
 */
Vec1D.prototype.c = function() {
	window.console.log(this);
	return this;
};

/**
 * Object for working with 2D vectors
 * @param {number,Vec1D} x number or another Vec2D; if empty then assumed 0
 * @param {number} y       number, if empty then assumed 0
 */
const Vec2D = function(x, y) {
	if (x instanceof Vec2D) {
		this.x = x.x;
		this.y = x.y;
	} else {
		// x
		if (typeof x !== "undefined") {
			if (typeof x === "number") {
				this.x = x;
			} else {
				throw new TypeError("Vec2D, constructor: Invalid parameter 'x': must be a Vec2D or a number");
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof y !== "undefined") {
			if (typeof y === "number") {
				this.y = y;
			} else {
				throw new TypeError("Vec2D, constructor: Invalid parameter 'y': must be a number");
			}
		} else {
			this.y = 0.0;
		}
	}
};

/**
 * Add another vector to this vector
 * @param {Vec2D} v vector Vec2D to be added
 * @return {Vec2D}  new Vec2D instance
 * @throws {TypeError} If v is not a Vec2D
 */
Vec2D.prototype.add = function(v) {
	if (v instanceof Vec2D) {
		return new Vec2D(this.x + v.x, this.y + v.y);
	} else {
		throw new TypeError("Vec2D, add: Invalid parameter: must be a Vec2D");
	}
};

/**
 * Multiplication by a scalar or by a vector
 * @param  {number,Vec2D} m a scalar or a Vec2D
 * @return {Vec2D}          new Vec2D instance
 * @throws {TypeError} If m is not a number or a Vec2D
 */
Vec2D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec2D(this.x * m, this.y * m);
	} else if (m instanceof Vec2D) {
		return new Vec2D(this.x * m.x, this.y * m.y);
	} else {
		throw new TypeError("Vec2D, mul: Invalid parameter: must be a number or a Vec2D");
	}
};

/**
 * Dot product of this Vec2D and another Vec2D
 * @param  {Vec2D} v another Vec2D
 * @return {number}  dot product
 * @throws {TypeError} If v is not a Vec2D
 */
Vec2D.prototype.dot = function(v) {
	if (v instanceof Vec2D) {
		return this.x * v.x + this.y * v.y;
	} else {
		throw new TypeError("Vec2D, dot: Invalid parameter: must be a Vec2D");
	}
};

/**
 * Length of this Vec2D vector
 * @return {number} length
 */
Vec2D.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y);
};

/**
 * Vector normalization
 * @return {Vec2D} new instance of a normalized Vec2D
 */
Vec2D.prototype.normalized = function() {
	const len = this.length();
	if (len === 0.0) return new Vec2D(0, 0);
	return new Vec2D(this.x / len, this.y / len);
};

/**
 * Print values to the console
 * @return {Vec2D} reference to this instance for calls chaining
 */
Vec2D.prototype.c = function() {
	window.console.log(this);
	return this;
};


/**
 * Object for working with 3D vectors
 * @param {number,Vec3D,Point3D} x a number or a Vec3D or a Point3D, if empty then assumed 0
 * @param {number} y               number, if empty then assumed 0
 * @param {number} z               number, if empty then assumed 0
 * @constructor
 */
const Vec3D = function(x, y, z) {
	if (x instanceof Vec3D || x instanceof Point3D) {
		this.x = x.x;
		this.y = x.y;
		this.z = x.z;
	} else {
		// x
		if (typeof x !== "undefined") {
			if (typeof x === "number") {
				this.x = x;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'x': must be a Vec3D, a Point3D or a number");
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof y !== "undefined") {
			if (typeof y === "number") {
				this.y = y;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'y': must be a number");
			}
		} else {
			this.y = 0.0;
		}
		// z
		if (typeof z !== "undefined") {
			if (typeof z === "number") {
				this.z = z;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'z': must be a number");
			}
		} else {
			this.z = 0.0;
		}
	}
};

/**
 * Add another vector to this vector
 * @param {Vec3D} v vector Vec3D to be added
 * @return {Vec3D}  new Vec3D instance
 * @throws {TypeError} If v is not a Vec3D
 */
Vec3D.prototype.add = function(v) {
	if (v instanceof Vec3D) {
		return new Vec3D(this.x + v.x, this.y + v.y, this.z + v.z);
	} else {
		throw new TypeError("Vec3D, add: Invalid parameter: must be a Vec3D");
	}
};

/**
 * Multiplication by a scalar or by a vector or by a quaternion or by a matrix (3×3) from right
 * @param  {number,Mat3,Vec3D,Quat} m a scalar or a Mat3 or a Quat or a Vec3D
 * @return {Vec3D}                    new Vec3D instance
 * @throws {TypeError} If m is not of a valid type
 */
Vec3D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec3D(this.x * m, this.y * m, this.z * m);
	} else if (m instanceof Mat3) {
		const temp = new Vec3D();
		temp.x = m.mat[0][0] * this.x + m.mat[1][0] * this.y + m.mat[2][0] * this.z;
		temp.y = m.mat[0][1] * this.x + m.mat[1][1] * this.y + m.mat[2][1] * this.z;
		temp.z = m.mat[0][2] * this.x + m.mat[1][2] * this.y + m.mat[2][2] * this.z;
		return temp;
	} else if (m instanceof Vec3D) {
		return new Vec3D(this.x * m.x, this.y * m.y, this.z * m.z);
	} else if (m instanceof Quat) {
		let p = new Quat(0, this.x, this.y, this.z);
		p = m.mulR(p).mulR(m.inv());
		return new Vec3D(p.i, p.j, p.k);
	} else {
		throw new TypeError("Vec3D, mul: Invalid parameter: must be a number or a Mat3 or a Quat or a Vec3D");
	}
};

/**
 * Dot product of this Vec3D and another Vec3D
 * @param  {Vec3D} v another Vec3D
 * @return {number}  dot product
 * @throws {TypeError} If v is not a Vec3D
 */
Vec3D.prototype.dot = function(v) {
	if (v instanceof Vec3D) {
		return this.x * v.x + this.y * v.y + this.z * v.z;
	} else {
		throw new TypeError("Vec3D, dot: Invalid parameter: must be a Vec3D");
	}
};

/**
 * Cross product of this Vec3D and another Vec3D
 * @param  {Vec3D} v another Vec3D
 * @return {Vec3D}   new Vec3D instance
 * @throws {TypeError} If v is not a Vec3D
 */
Vec3D.prototype.cross = function(v) {
	if (v instanceof Vec3D) {
		return new Vec3D(this.y * v.z - this.z * v.y, this.z * v.x - this.x * v.z, this.x * v.y - this.y * v.x);
	} else {
		throw new TypeError("Vec3D, cross: Invalid parameter: must be a Vec3D");
	}
};

/**
 * Length of this Vec3D vector
 * @return {number} length
 */
Vec3D.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
};

/**
 * Vector normalization
 * @return {Vec3D} new instance of a normalized Vec3D
 */
Vec3D.prototype.normalized = function() {
	const len = this.length();
	if (len === 0) return new Vec3D(0, 0, 0);
	return new Vec3D(this.x / len, this.y / len, this.z / len);
};

/**
 * Print values to the console
 * @return {Vec3D} reference to this instance for calls chaining
 */
Vec3D.prototype.c = function() {
	window.console.log(this);
	return this;
};


/**
 * Object for working with 3D points (points with homogeneous coordinate)
 * @param {number,Vec3D,Point3D} ax a number or a Vec3D or a Point3D, if empty then assumed 0
 * @param {number} ay               number, if empty then assumed 0
 * @param {number} az               number, if empty then assumed 0
 * @param {number} aw               number, if empty then assumed 1
 * @constructor
 */
const Point3D = function(ax, ay, az, aw) {
	if (ax instanceof Vec3D) {
		this.x = ax.x;
		this.y = ax.y;
		this.z = ax.z;
		this.w = 1.0;
	} else if (ax instanceof Point3D) {
		this.x = ax.x;
		this.y = ax.y;
		this.z = ax.z;
		this.w = ax.w;
	} else {
		// x
		if (typeof ax !== "undefined") {
			if (typeof ax === "number") {
				this.x = ax;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'ax': must be a Vec3D or a Point3D or a number");
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof ay !== "undefined") {
			if (typeof ay === "number") {
				this.y = ay;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'ay': must be a number");
			}
		} else {
			this.y = 0.0;
		}
		// z
		if (typeof az !== "undefined") {
			if (typeof az === "number") {
				this.z = az;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'az': must be a number");
			}
		} else {
			this.z = 0.0;
		}
		// w
		if (typeof aw !== "undefined") {
			if (typeof aw === "number") {
				this.w = aw;
			} else {
				throw new TypeError("Vec3D, constructor: Invalid parameter 'aw': must be a number");
			}
		} else {
			this.w = 1.0;
		}
	}
};

/**
 * Multiplication by a scalar or by a quaternion or by a matrix (4×4) from right
 * @param  {number,Mat4,Quat} m a scalar or a Mat4 or a Quat
 * @return {Point3D}            new Point3D instance
 * @throws {TypeError} If m is not of a valid type
 */
Point3D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Point3D(this.x * m, this.y * m, this.z * m, this.w * m);
	} else if (m instanceof Mat4) {
		const res = new Point3D();
		res.x = m.mat[0][0] * this.x + m.mat[1][0] * this.y + m.mat[2][0] * this.z + m.mat[3][0] * this.w;
		res.y = m.mat[0][1] * this.x + m.mat[1][1] * this.y + m.mat[2][1] * this.z + m.mat[3][1] * this.w;
		res.z = m.mat[0][2] * this.x + m.mat[1][2] * this.y + m.mat[2][2] * this.z + m.mat[3][2] * this.w;
		res.w = m.mat[0][3] * this.x + m.mat[1][3] * this.y + m.mat[2][3] * this.z + m.mat[3][3] * this.w;
		return res;
	} else if (m instanceof Quat) {
		return new Point3D(this.dehomog().mul(m));
	} else {
		throw new TypeError("Point3D, mul: Invalid parameter: must be a number or a Mat4 or a Quat");
	}
};

/**
 * Add another point to this point
 * @param {Point3D} p point Point3D to be added
 * @return {Point3D}  new Point3D instance
 * @throws {TypeError} If v is not a Point3D
 */
Point3D.prototype.add = function(p) {
	if (p instanceof Point3D) {
		return new Point3D(this.x + p.x, this.y + p.y, this.z + p.z, this.w + p.w);
	} else {
		throw new TypeError("Point3D, add: Invalid parameter: must be a Point3D");
	}
};

/**
 * Subtract another point from this point
 * @param {Point3D} p point Point3D to be subtracted
 * @return {Point3D}  new Point3D instance
 * @throws {TypeError} If v is not a Point3D
 */
Point3D.prototype.sub = function(p) {
	if (p instanceof Point3D) {
		return new Point3D(this.x - p.x, this.y - p.y, this.z - p.z, this.w - p.w);
	} else {
		throw new TypeError("Point3D, sub: Invalid parameter: must be a Point3D");
	}
};

/**
 * Dehomogenization of this point (divide by w coordinate)
 * @return {Vec3D} new Vec3D instance
 */
Point3D.prototype.dehomog = function() {
	if (this.w === 0) return new Vec3D(0, 0, 0);
	return new Vec3D(this.x / this.w, this.y / this.w, this.z / this.w);
};

/**
 * Transform to Vec3D, ignore w coordinate
 * @return {Vec3D} new Vec3D instance
 */
Point3D.prototype.ignoreW = function() {
	return new Vec3D(this.x, this.y, this.z);
};

/**
 * Print values to the console
 * @return {Point3D} reference to this instance for calls chaining
 */
Point3D.prototype.c = function() {
	window.console.log(this);
	return this;
};

/**
 * Object for working with quaternions
 * @param {number,Vec3D,Point3D} r a number or a Quat, if empty then assumed 0
 * @param {number} i               number, if empty then assumed 0
 * @param {number} j               number, if empty then assumed 0
 * @param {number} k               number, if empty then assumed 0
 * @constructor
 */
const Quat = function(r, i, j, k) {
	if (r instanceof Quat) {
		this.i = r.i;
		this.j = r.j;
		this.k = r.k;
		this.r = r.r;
	} else {
		// i
		if (typeof i !== "undefined") {
			if (typeof i === "number") {
				this.i = i;
			} else {
				throw new TypeError("Quat, constructor: Invalid parameter i: must be a number or a Quat");
			}
		} else {
			this.i = 0.0;
		}
		// j
		if (typeof j !== "undefined") {
			if (typeof j === "number") {
				this.j = j;
			} else {
				throw new TypeError("Quat, constructor: Invalid parameter j: must be a number");
			}
		} else {
			this.j = 0.0;
		}
		// k
		if (typeof k !== "undefined") {
			if (typeof k === "number") {
				this.k = k;
			} else {
				throw new TypeError("Quat, constructor: Invalid parameter k: must be a number");
			}
		} else {
			this.k = 0.0;
		}
		// r
		if (typeof r !== "undefined") {
			if (typeof r === "number") {
				this.r = r;
			} else {
				throw new TypeError("Quat, constructor: Invalid parameter r: must be a number");
			}
		} else {
			this.r = 0.0;
		}
	}
};

/**
 * Add another quaternion to this quaternion
 * @param {Quat} q quaternion to be added
 * @return {Quat}  new Quat instance
 * @throws {TypeError} If q is not a Quat
 */
Quat.prototype.add = function(q) {
	if (q instanceof Quat) {
		return new Quat(this.r + q.r, this.i + q.i, this.j + q.j, this.k + q.k);
	} else {
		throw new TypeError("Quat, add: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns the result of quaternion subtraction of the given quaternion
 * @param {Quat} q quaternion to subtract
 * @return {Quat}  new Quat instance
 * @throws {TypeError} If q is not a Quat
 */
Quat.prototype.sub = function(q) {
	if (q instanceof Quat) {
		return new Quat(this.r - q.r, this.i - q.i, this.j - q.j, this.k - q.k);
	} else {
		throw new TypeError("Quat, sub: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns the result of right side quaternion multiplication by the given quaternion
 * @param {Quat} q quaternion to multiplicate
 * @return {Quat}  new Quat instance
 * @throws {TypeError} If q is not a Quat
 */
Quat.prototype.mulR = function(q) {
	if (q instanceof Quat) {
		return new Quat(
			this.r * q.r - this.i * q.i - this.j * q.j - this.k * q.k,
			this.r * q.i + this.i * q.r + this.j * q.k - this.k * q.j,
			this.r * q.j - this.i * q.k + this.j * q.r + this.k * q.i,
			this.r * q.k + this.i * q.j - this.j * q.i + this.k * q.r
		);
	} else {
		throw new TypeError("Quat, mulR: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns the result of left side quaternion multiplication by the given quaternion
 * @param {Quat} q quaternion to be multiplicated
 * @return {Quat}  new Quat instance
 * @throws {TypeError} If q is not a Quat
 */
Quat.prototype.mulL = function(q) {
	if (q instanceof Quat) {
		return new Quat(
			q.r * this.r - q.i * this.i - q.j * this.j - q.k * this.k,
			q.r * this.i + q.i * this.r + q.j * this.k - q.k * this.j,
			q.r * this.j + q.j * this.r + q.k * this.i - q.i * this.k,
			q.r * this.k + q.k * this.r + q.i * this.j - q.j * this.i
		);
	} else {
		throw new TypeError("Quat, mulL: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns the result of scalar multiplication or right side quaternion multiplication
 * @param {number,Quat} q a scalar or a Quat
 * @return {Quat}         new Quat instance
 * @throws {TypeError} If q is not a number or a Quat
 */
Quat.prototype.mul = function(q) {
	if (typeof q === "number") {
		return new Quat(q * this.r, q * this.i, q * this.j, q * this.k);
	} else if (q instanceof Quat) {
		return this.mulR(q);
	} else {
		throw new TypeError("Quat, mul: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns the norm of this quaternion
 * @return {number} norm
 */
Quat.prototype.norma = function() {
	return Math.sqrt(this.r * this.r + this.i * this.i + this.j * this.j + this.k * this.k);
};

/**
 * Returns the inverse of this quaternion if it exists or an empty quaternion
 * @return {Quat} new Quat instance
 */
Quat.prototype.inv = function() {
	const n = this.norma();
	const norm = n * n;
	if (norm > 0) {
		return new Quat(this.r / norm, -this.i / norm, -this.j / norm, -this.k / norm);
	} else {
		return new Quat(0, 0, 0, 0);
	}
};

/**
 * Returns logarithm function of this quaternion
 * @return {Quat} new Quat instance
 */
Quat.prototype.log = function() {
	if (this.i === 0 && this.j === 0 && this.k === 0) {
		if (this.r > 0) {
			return new Quat(Math.log(this.r), 0, 0, 0);
		} else if (this.r < 0) {
			return new Quat(Math.log(-this.r), 1, 0, 0);
		} else {
			return new Quat();
		}
	} else {
		const s = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
		const a = Math.atan2(s, this.r) / s;
		return new Quat(Math.log(this.norma()), a * this.i, a * this.j, a * this.k);
	}
};

/**
 * Returns exponential function of this quaternion
 * @return {Quat} new Quat instance
 */
Quat.prototype.exp = function() {
	if (this.i === 0 && this.j === 0 && this.k === 0) {
		return new Quat(Math.exp(this.r), 0, 0, 0);
	} else {
		const s1 = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
		const cos = Math.cos(s1);
		const s = Math.exp(this.r) * Math.sin(s1) / s1;
		return new Quat(Math.exp(this.r) * cos, s * this.i, s * this.j, s * this.k);
	}
};

/**
 * Returns the opposite quaternion to this quaternion
 * @return {Quat} new Quat instance
 */
Quat.prototype.neg = function() {
	return new Quat(-this.r, -this.i, -this.j, -this.k);
};

/**
 * Returns the result of dot-product with the given quaternion
 * @param  {Quat} q  a quaternion
 * @return {number}  a number
 */
Quat.prototype.dot = function(q) {
	if (q instanceof Quat) {
		return this.i * q.i + this.j * q.j + this.k * q.k + this.r * q.r;
	} else {
		throw new TypeError("Quat, dot: Invalid parameter: must be a Quat");
	}
};

/**
 * Returns a normalized quaternion if possible (nonzero norm), empty quaternion otherwise
 * @return {Quat} new Quat instance
 */
Quat.prototype.renorm = function() {
	const norm = this.norma();
	if (norm > 0) {
		return new Quat(this.r / norm, this.i / norm, this.j / norm, this.k / norm);
	} else {
		return new Quat(0, 0, 0, 0);
	}
};

/**
 * Creates a 4×4 transformation matrix equivalent to rotation defined by this quaternion
 * @return {Mat4} new Mat4 instance
 */
Quat.prototype.toRotationMatrix = function() {
	const res = new Mat4Identity();
	this.renorm();
	res.mat[0][0] = 1 - 2 * (this.j * this.j + this.k * this.k);
	res.mat[1][0] = 	2 * (this.i * this.j - this.r * this.k);
	res.mat[2][0] = 	2 * (this.r * this.j + this.i * this.k);

	res.mat[0][1] = 	2 * (this.i * this.j + this.r * this.k);
	res.mat[1][1] = 1 - 2 * (this.i * this.i + this.k * this.k);
	res.mat[2][1] = 	2 * (this.k * this.j - this.i * this.r);

	res.mat[0][2] = 	2 * (this.i * this.k - this.r * this.j);
	res.mat[1][2] = 	2 * (this.k * this.j + this.i * this.r);
	res.mat[2][2] = 1 - 2 * (this.i * this.i + this.j * this.j);
	return res;
};

/**
 * Helper object for work with quaternions
 * @static
 */
const Quat2 = {};

/**
 * Creates a new quaternion equivalent to the rotation given by the 4×4 transformation matrix
 * @param  {Mat4} mat input rotation matrix
 * @return {Quat}     new Quat instance
 * @throws {TypeError} If mat is not a Mat4
 */
Quat2.fromRotationMatrix = function(mat) {
	if (mat instanceof Mat4) {
		let r, i, j, k;
		const diagonal1 = mat.mat[0][0] + mat.mat[1][1] + mat.mat[2][2];

		if (diagonal1 > 0.0) {
			r = (0.5 * Math.sqrt(diagonal1 + mat.mat[3][3]));
			i = (mat.mat[2][1] - mat.mat[1][2]) / (4 * r);
			j = (mat.mat[0][2] - mat.mat[2][0]) / (4 * r);
			k = (mat.mat[1][0] - mat.mat[0][1]) / (4 * r);
		} else {
			const indices = [1, 2, 0];
			let a = 0;

			if (mat.mat[1][1] > mat.mat[0][0]) a = 1;
			if (mat.mat[2][2] > mat.mat[a][a]) a = 2;

			const b = indices[a];
			const c = indices[b];

			const diagonal2 = mat.mat[a][a] - mat.mat[b][b] - mat.mat[c][c] + mat.mat[3][3];
			r = (0.5 * Math.sqrt(diagonal2));
			i = (mat.mat[a][b] + mat.mat[b][a]) / (4 * r);
			j = (mat.mat[a][c] + mat.mat[c][a]) / (4 * r);
			k = (mat.mat[b][c] - mat.mat[c][b]) / (4 * r);
		}
		return new Quat(r, i, j, k);
	} else {
		throw new TypeError("Quat, fromRotationMatrix: Invalid parameter: must be a Mat4");
	}
};

/**
 * Creates a new quaternion equivalent to the rotations around axis given by vector x, y and z
 * @param  {number} angle rotation angle in radians
 * @param  {number} a     x vector coordinate
 * @param  {number} b     y vector coordinate
 * @param  {number} c     z vector coordinate
 * @return {Quat}         new Quat instance
 * @throws {TypeError} If any of the parameters is not a number
 */
Quat2.fromEulerAngle = function(angle, a, b, c) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, fromEulerAngle: Invalid number of parameters: must be 4");
	} else if (typeof angle === "number" && typeof a === "number" && typeof b === "number" && typeof c === "number") {
		return new Quat(
			Math.cos(angle / 2),
			Math.sin(angle / 2) * a,
			Math.sin(angle / 2) * b,
			Math.sin(angle / 2) * c
		);
	} else {
		throw new TypeError("Quat, fromEulerAngle: Invalid parameter: must be 4 numbers");
	}
};

/**
 * Creates a new quaternion equivalent to the right-handed rotations
 *     given by the Euler angles around x, y and z axes chained in sequence in this order
 * @param  {number} a rotation angle around x-axis
 * @param  {number} b rotation angle around y-axis
 * @param  {number} c rotation angle around z-axis
 * @return {Quat}     new Quat instance
 * @throws {TypeError} If any of the parameters is not a number
 */
Quat2.fromEulerAngles = function(a, b, c) {
	if (arguments.length !== 3) {
		throw new TypeError("Quat, fromEulerAngles: Invalid number of parameters: must be 3");
	} else if (typeof a === "number" && typeof b === "number" && typeof c === "number") {
		const Qi = Quat2.fromEulerAngle(a, 1, 0, 0);
		const Qj = Quat2.fromEulerAngle(b, 0, 1, 0);
		const Qk = Quat2.fromEulerAngle(c, 0, 0, 1);
		return new Quat(Qk.mul(Qj).mul(Qi));
	} else {
		throw new TypeError("Quat, fromEulerAngles: Invalid parameter: must be 3 numbers");
	}
};

/**
 * Returns the angle and axis rotation in format Point3D:(angle, x-axis, y-axis, z-axis)
 * @return {Point3D} new Point3D instance
 */
Quat.prototype.toEulerAngle = function() {
	const angle = 2 * Math.acos(this.r);
	const x = this.i;
	const y = this.j;
	const z = this.k;

	const s = Math.sqrt(x * x + y * y + z * z);
	if (s < 0.0001) {
		return new Point3D(angle, 1.0, 0.0, 0.0);
	} else {
		return new Point3D(angle, (x / s), (y / s),	(z / s));
	}
};

/**
 * Linear interpolation between this and another quaternion Lerp(Q1,Q2,t)=(1-t)Q1+tQ2
 * @param  {Quat} q   another quaternion
 * @param  {number} t interpolation parameter in interval <0;1>
 * @return {Quat}     new Quat instance
 * @throws {TypeError} If q is not a Quat and t is not a number
 */
Quat.prototype.lerp = function(q, t) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, lerp: Invalid number of parameters: must be 2");
	} else if (q instanceof Quat && typeof t === "number") {
		if (t >= 1) {
			return new Quat(q);
		} else if (t <= 0) {
			return new Quat(this);
		} else {
			return new Quat((this.mul(1 - t)).add(q.mul(t)));
		}
	} else {
		throw new TypeError("Quat, lerp: Invalid parameter: must be a Quat and a number");
	}
};

/**
 * Spherical interpolation between this and another quaternion
 * @param  {Quat} q   another quaternion
 * @param  {number} t interpolation parameter in interval <0;1>
 * @return {Point3D} new Point3D instance
 * @throws {TypeError} If q is not a Quat and t is not a number
 */
Quat.prototype.slerp = function(q, t) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, lerp: Invalid number of parameters: must be 2");
	} else if (q instanceof Quat && typeof t === "number") {
		let c = this.dot(q);
		if (c > 1.0) {
			c = 1.0;
		} else if (c < -1.0) {
			c = -1.0;
		}
		const angle = Math.acos(c);
		if (Math.abs(angle) < 1.0e-5) {
			return new Quat(this);
		}
		const s = 1 / Math.sin(angle);
		if (t >= 1) {
			return new Quat(this);
		} else if (t <= 0) {
			return new Quat(q);
		} else {
			return new Quat(
				this.renorm()
					.mul(Math.sin((1 - t) * angle) * s)
					.add(q.renorm().mul(Math.sin(t * angle) * s))
				).renorm();
		}
	} else {
		throw new TypeError("Quat, lerp: Invalid parameter: must be a Quat and a number");
	}
};

/**
 * Cubic interpolation between this and another quaternion
 * @param  {Quat} q   another quaternion
 * @param  {Quat} q1  another quaternion
 * @param  {Quat} q2  another quaternion
 * @param  {number} t interpolation parameter in interval <0;1>
 * @return {Quat}     new Quat instance
 * @throws {TypeError} If q, q1, q2 are not Quats and t is not a number
 */
Quat.prototype.squad = function(q, q1, q2, t) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, squad: Invalid number of parameters: must be 4");
	} else if (q instanceof Quat && q1 instanceof Quat && q2 instanceof Quat && typeof t === "number") {
		return new Quat(this.slerp(q, t).slerp(q1.slerp(q2, t), (2 * t * (1 - t))));
	} else {
		throw new TypeError("Quat, squad: Invalid parameter: must be 3 Quats and a number");
	}
};

/**
 *
 * @param  {Quat} q1 quaternion
 * @param  {Quat} q2 quaternion
 * @return {Quat}    new Quat instance
 * @throws {TypeError} If q1 and q2 are not Quats
 */
Quat.prototype.quadrangle = function(q1, q2) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, quadrangle: Invalid number of parameters: must be 2");
	} else if (q1 instanceof Quat && q2 instanceof Quat) {
		const s1 = this.inv().mul(q1);
		const s2 = this.inv().mul(q2);
		return new Quat((s1.log().add(s2.log()).mul(-1 / 4)).exp());
	} else {
		throw new TypeError("Quat, quadrangle: Invalid parameter: must be 2 Quats");
	}
};

/**
 *
 * @param  {Quat} q1   quaternion
 * @param  {Quat} q2   quaternion
 * @param  {Quat} q3   quaternion
 * @param  {number}  t interpolation parameter in interval <0;1>
 * @return {Quat}      new Quat instance
 * @throws {TypeError} If q1, q2, q3 are not Quats and t is not a number
 */
Quat.prototype.squad2 = function(q1, q2, q3, t) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, squad2: Invalid number of parameters: must be 4");
	} else if (q1 instanceof Quat && q2 instanceof Quat && q3 instanceof Quat && typeof t === "number") {
		const s1 = this.quadrangle(q1, q2);
		const s2 = q2.quadrangle(this, q3);
		return new Quat(this.slerp(q2, t).slerp(s1.slerp(s2, t), (2 * t * (1 - t))));
	} else {
		throw new TypeError("Quat, squad2: Invalid parameter: must be 3 Quats and a number");
	}
};

/**
 * Print values to the console
 * @return {Quat} reference to this instance for calls chaining
 */
Quat.prototype.c = function() {
	window.console.log(this);
	return this;
};


/**
 * Object for work with 3×3 matrices
 * @param {Vec3D,Mat3,Mat4} v1 Vec3D or Mat3 or Mat4
 * @param {Vec3D} v2           Vec3D
 * @param {Vec3D} v3           Vec3D
 * @constructor
 */
const Mat3 = function(v1, v2, v3) {
	this.mat = [];
	if (v1 instanceof Vec3D && v2 instanceof Vec3D && v3 instanceof Vec3D) {
		this.mat[0] = [];
		this.mat[0][0] = v1.x;
		this.mat[0][1] = v1.y;
		this.mat[0][2] = v1.z;
		this.mat[1] = [];
		this.mat[1][0] = v2.x;
		this.mat[1][1] = v2.y;
		this.mat[1][2] = v2.z;
		this.mat[2] = [];
		this.mat[2][0] = v3.x;
		this.mat[2][1] = v3.y;
		this.mat[2][2] = v3.z;
	} else if (v1 instanceof Mat3) {
		for (let i = 0; i < 3; i++) {
			this.mat[i] = [];
			for (let j = 0; j < 3; j++) {
				this.mat[i][j] = v1.mat[i][j];
			}
		}
	} else if (v1 instanceof Mat4) {
		for (let i = 0; i < 3; i++) {
			this.mat[i] = [];
			for (let j = 0; j < 3; j++) {
				this.mat[i][j] = v1.mat[i][j];
			}
		}
	} else {
		this.mat = new ZeroArray(3);
	}
};

/**
 * Creates 3×3 identity matrix
 * @augments {Mat3}
 */
const Mat3Identity = function() {
	this.mat = new ZeroArray(3);
	for (let i = 0; i < 3; i++) {
		this.mat[i][i] = 1;
	}
};
Mat3Identity.prototype = Object.create(Mat3.prototype);
Mat3Identity.prototype.constructor = Mat3Identity;
Mat3Identity.prototype.parent = Mat3;

/**
 * Creates a transformation matrix 3×3 to rotate around the X axis in 3D
 * @augments {Mat3}
 * @param {number} alpha rotation angle in radians
 */
const Mat3RotX = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotX, constructor: Invalid parameter: must be a number");
	}
	this.mat = new Mat3Identity().mat;
	this.mat[1][1] = Math.cos(alpha);
	this.mat[2][2] = Math.cos(alpha);
	this.mat[2][1] = -Math.sin(alpha);
	this.mat[1][2] = Math.sin(alpha);
};
Mat3RotX.prototype = Object.create(Mat3.prototype);
Mat3RotX.prototype.constructor = Mat3RotX;
Mat3RotX.prototype.parent = Mat3;

/**
 * Creates a transformation matrix 3×3 to rotate around the Y axis in 3D
 * @augments {Mat3}
 * @param {number} alpha rotation angle in radians
 */
const Mat3RotY = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotY, constructor: Invalid parameter: must be a number");
	}
	this.mat = new Mat3Identity().mat;
	this.mat[0][0] = Math.cos(alpha);
	this.mat[2][2] = Math.cos(alpha);
	this.mat[2][0] = Math.sin(alpha);
	this.mat[0][2]= -Math.sin(alpha);
};
Mat3RotY.prototype = Object.create(Mat3.prototype);
Mat3RotY.prototype.constructor = Mat3RotY;
Mat3RotY.prototype.parent = Mat3;

/**
 * Creates a transformation matrix 3×3 to rotate around the Z axis in 3D
 * @augments {Mat3}
 * @param {number} alpha rotation angle in radians
 */
const Mat3RotZ = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotZ, constructor: Invalid parameter: must be a number");
	}
	this.mat = new Mat3Identity().mat;
	this.mat[0][0] = Math.cos(alpha);
	this.mat[1][1] = Math.cos(alpha);
	this.mat[1][0] = -Math.sin(alpha);
	this.mat[0][1] = Math.sin(alpha);
};
Mat3RotZ.prototype = Object.create(Mat3.prototype);
Mat3RotZ.prototype.constructor = Mat3RotZ;
Mat3RotZ.prototype.parent = Mat3;

/**
 * Returns the result of element-wise summation with another 3×3 matrix
 * @param {Mat3} m 3×3 matrix
 * @return {Mat3}  new Mat3 instance
 * @throws {TypeError} If m is not a Mat3
 */
Mat3.prototype.add = function(m) {
	if (m instanceof Mat3) {
		const temp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				temp.mat[i][j] = this.mat[i][j] + m.mat[i][j];
			}
		}
		return temp;
	} else {
		throw new TypeError("Mat3, add: Invalid parameter: must be a Mat3");
	}
};

/**
 * Returns the result of element-wise multiplication by the given scalar value or the result of matrix multiplication
 * @param  {number,Mat3} m a number or a Mat3
 * @return {Mat3}          new Mat3 instance
 * @throws {TypeError} If m is not a number or a Mat3
 */
Mat3.prototype.mul = function(m) {
	if (typeof m === "number") {
		const temp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				temp.mat[i][j] = this.mat[i][j] * m;
			}
		}
		return temp;
	} else if (m instanceof Mat3) {
		const temp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				let sum = 0;
				for (let k = 0; k < 3; k++) {
					sum += this.mat[i][k] * m.mat[k][j];
				}
				temp.mat[i][j] = sum;
			}
		}
		return temp;
	} else {
		throw new TypeError("Mat3, mul: Invalid parameter: must be a number nebo a Mat3");
	}
};

/**
 * Returns the transposition of this matrix
 * @return {Mat3} new Mat3 instance
 */
Mat3.prototype.transpose = function() {
	const temp = new Mat3();
	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			temp.mat[i][j] = this.mat[j][i];
		}
	}
	return temp;
};

/**
 * Returns the determinant of this matrix
 * @return {number} determinant
 */
Mat3.prototype.det = function() {
	return this.mat[0][0] * (this.mat[1][1] * this.mat[2][2] - this.mat[2][1] * this.mat[1][2]) -
	       this.mat[0][1] * (this.mat[1][0] * this.mat[2][2] - this.mat[2][0] * this.mat[1][2]) +
	       this.mat[0][2] * (this.mat[1][0] * this.mat[2][1] - this.mat[2][0] * this.mat[1][1]);
};

/**
 * Returns the inverse of this matrix if it exists
 * @return {Mat3} new Mat3 instance
 */
Mat3.prototype.inverse = function() {
	const temp = new Mat3();
	const det = 1.0 / this.det();
	temp.mat[0][0] = det * (this.mat[1][1] * this.mat[2][2] - this.mat[1][2] * this.mat[2][1]);
	temp.mat[0][1] = det * (this.mat[0][2] * this.mat[2][1] - this.mat[0][1] * this.mat[2][2]);
	temp.mat[0][2] = det * (this.mat[0][1] * this.mat[1][2] - this.mat[0][2] * this.mat[1][1]);

	temp.mat[1][0] = det * (this.mat[1][2] * this.mat[2][0] - this.mat[1][0] * this.mat[2][2]);
	temp.mat[1][1] = det * (this.mat[0][0] * this.mat[2][2] - this.mat[0][2] * this.mat[2][0]);
	temp.mat[1][2] = det * (this.mat[0][2] * this.mat[1][0] - this.mat[0][0] * this.mat[1][2]);

	temp.mat[2][0] = det * (this.mat[1][0] * this.mat[2][1] - this.mat[1][1] * this.mat[2][0]);
	temp.mat[2][1] = det * (this.mat[0][1] * this.mat[2][0] - this.mat[0][0] * this.mat[2][1]);
	temp.mat[2][2] = det * (this.mat[0][0] * this.mat[1][1] - this.mat[0][1] * this.mat[1][0]);
	return temp;
};

/**
 * Print values to the console
 * @return {Mat3} reference to this instance for calls chaining
 */
Mat3.prototype.c = function() {
	window.console.log("%cMat3", 'font-style: italic; font-weight: bold');
	for (let i = 0; i < 3; i++) {
		let x = "";
		for (let j = 0; j < 3; j++) {
			x += this.mat[i][j] + ", ";
		}
		window.console.log(x.substring(0, x.length - 2));
	}
	return this;
};


/**
 * Object for work with 4×4 matrices
 * @param {Point3D,Mat4} p1 Point3D or Mat4
 * @param {Point3D} p2      Point3D
 * @param {Point3D} p3      Point3D
 * @param {Point3D} p4      Point3D
 * @constructor
 */
const Mat4 = function(p1, p2, p3, p4) {
	this.mat = [];
	if (p1 instanceof Point3D && p2 instanceof Point3D && p3 instanceof Point3D && p4 instanceof Point3D) {
		this.mat[0] = [];
		this.mat[0][0] = p1.x;
		this.mat[0][1] = p1.y;
		this.mat[0][2] = p1.z;
		this.mat[0][3] = p1.w;
		this.mat[1] = [];
		this.mat[1][0] = p2.x;
		this.mat[1][1] = p2.y;
		this.mat[1][2] = p2.z;
		this.mat[1][3] = p2.w;
		this.mat[2] = [];
		this.mat[2][0] = p3.x;
		this.mat[2][1] = p3.y;
		this.mat[2][2] = p3.z;
		this.mat[2][3] = p3.w;
		this.mat[3] = [];
		this.mat[3][0] = p4.x;
		this.mat[3][1] = p4.y;
		this.mat[3][2] = p4.z;
		this.mat[3][3] = p4.w;
	} else if (p1 instanceof Mat4) {
		for (let i = 0; i < 4; i++) {
			this.mat[i] = [];
			for (let j = 0; j < 4; j++) {
				this.mat[i][j] = p1.mat[i][j];
			}
		}
	} else {
		this.mat = new ZeroArray(4);
	}
};

/**
 * Creates 4×4 identity matrix
 * @augments {Mat4}
 */
const Mat4Identity = function() {
	this.mat = new ZeroArray(4);
	for (let i = 0; i < 4; i++) {
		this.mat[i][i] = 1.0;
	}
};
Mat4Identity.prototype = Object.create(Mat4.prototype);
Mat4Identity.prototype.constructor = Mat4Identity;
Mat4Identity.prototype.parent = Mat4;

/**
 * Creates matrix of orthogonal visibility volume to normalized clipping volume transformation
 * @augments {Mat4}
 * @param {number} w  visibility width (usually canvas width)
 * @param {number} h  visibility height (usually canvas height)
 * @param {number} zn z-near, distance to the near clipping plane along z-axis
 * @param {number} zf z-far,distance to the far clipping plane along z-axis
 * @throws {TypeError} If any of the parameters is not defined or is not a number
 */
const Mat4OrthoRH = function(w, h, zn, zf) {
	if (arguments.length !== 4) {
		throw new TypeError("Mat4OrthoRH: Invalid number of parameters: must be 4");
	}
	if (typeof w !== "number" || typeof h !== "number" || typeof zn !== "number" || typeof zf !== "number") {
		throw new TypeError("Mat4OrthoRH: Invalid parameter: must be 4 numbers");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[0][0] = 2.0 / w;
	this.mat[1][1] = 2.0 / h;
	this.mat[2][2] = 1.0 / (zn - zf);
	this.mat[3][0]= -1;
	this.mat[3][1] = -1;
	this.mat[3][2] = zn / (zn - zf);
};
Mat4OrthoRH.prototype = Object.create(Mat4.prototype);
Mat4OrthoRH.prototype.constructor = Mat4OrthoRH;
Mat4OrthoRH.prototype.parent = Mat4;

/**
 * Creates matrix of perspective visibility volume to normalized clipping volume transformation
 * @augments {Mat4}
 * @param {number} alpha vertical field of view angle in radians
 * @param {number} k     volume height/width ratio
 * @param {number} zn    z-near, distance to the near clipping plane along z-axis
 * @param {number} zf    z-far, distance to the far clipping plane along z-axis
 * @throws {TypeError} If any of the parameters is not defined or is not a number
 */
const Mat4PerspRH = function(alpha, k, zn, zf) {
	if (arguments.length !== 4) {
		throw new TypeError("Mat4PerspRH: Invalid number of parameters: must be 4");
	}
	if (typeof alpha !== "number" || typeof k !== "number" || typeof zn !== "number" || typeof zf !== "number") {
		throw new TypeError("Mat4PerspRH: Invalid parameter: must be 4 numbers");
	}

	const h = (1.0 / Math.tan(alpha / 2.0));
	const w = k * h;
	this.mat = new Mat4Identity().mat;
	this.mat[0][0] = w;
	this.mat[1][1] = h;
	this.mat[2][2] = zf / (zn - zf);
	this.mat[3][2] = zn * zf / (zn - zf);
	this.mat[2][3] = -1.0;

};
Mat4PerspRH.prototype = Object.create(Mat4.prototype);
Mat4PerspRH.prototype.constructor = Mat4PerspRH;
Mat4PerspRH.prototype.parent = Mat4;


/**
 * Creates a transformation matrix 4×4 to rotate around the X axis in 3D
 * @augments {Mat4}
 * @param {number} alpha rotation angle in radians
 * @throws {TypeError} If alpha is not a number
 */
const Mat4RotX = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotX: Invalid parameter: must be a number");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[1][1] = Math.cos(alpha);
	this.mat[2][2] = Math.cos(alpha);
	this.mat[2][1] = -Math.sin(alpha);
	this.mat[1][2] = Math.sin(alpha);
};
Mat4RotX.prototype = Object.create(Mat4.prototype);
Mat4RotX.prototype.constructor = Mat4RotX;
Mat4RotX.prototype.parent = Mat4;

/**
 * Creates a transformation matrix 4×4 to rotate around the Y axis in 3D
 * @augments {Mat4}
 * @param {number} alpha rotation angle in radians
 * @throws {TypeError} If alpha is not a number
 */
const Mat4RotY = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotY: Invalid parameter: must be a number");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[0][0] = Math.cos(alpha);
	this.mat[2][2] = Math.cos(alpha);
	this.mat[2][0] = Math.sin(alpha);
	this.mat[0][2] = -Math.sin(alpha);
};
Mat4RotY.prototype = Object.create(Mat4.prototype);
Mat4RotY.prototype.constructor = Mat4RotY;
Mat4RotY.prototype.parent = Mat4;

/**
 * Creates a transformation matrix 4×4 to rotate around the Z axis in 3D
 * @augments {Mat4}
 * @param {number} alpha rotation angle in radians
 * @throws {TypeError} If alpha is not a number
 */
const Mat4RotZ = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotZ: Invalid parameter: must be a number");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[0][0] = Math.cos(alpha);
	this.mat[1][1] = Math.cos(alpha);
	this.mat[1][0] = -Math.sin(alpha);
	this.mat[0][1] = Math.sin(alpha);
};
Mat4RotZ.prototype = Object.create(Mat4.prototype);
Mat4RotZ.prototype.constructor = Mat4RotZ;
Mat4RotZ.prototype.parent = Mat4;

/**
 * Creates a transformation matrix 4×4 to rotate around the X,Y,Z axis in 3D
 * @augments {Mat4}
 * @param {number} alpha rotation angle in radians around X-axis
 * @param {number} beta  rotation angle in radians around Y-axis
 * @param {number} gamma rotation angle in radians around Z-axis
 * @throws {TypeError} If any of the parameters is not defined or is not a number
 */
const Mat4RotXYZ = function(alpha, beta, gamma) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4RotXYZ: Invalid number of parameters: must be 3");
	}
	if (typeof alpha !== "number" || typeof beta !== "number" || typeof gamma !== "number") {
		throw new TypeError("Mat4RotXYZ: Invalid parameter: must be 3 numbers");
	}

	this.mat = new Mat4RotX(alpha).mul(new Mat4RotY(beta)).mul(new Mat4RotZ(gamma)).mat;
};
Mat4RotXYZ.prototype = Object.create(Mat4.prototype);
Mat4RotXYZ.prototype.constructor = Mat4RotXYZ;
Mat4RotXYZ.prototype.parent = Mat4;

/**
 * Creates a transformation matrix 4×4 to scale in 3D
 * @augments {Mat4}
 * @param {number} x X-axis scale factor
 * @param {number} y Y-axis scale factor
 * @param {number} z Z-axis scale factor
 * @throws {TypeError} If any of the parameters is not defined or is not a number
 */
const Mat4Scale = function(x, y, z) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4Scale: Invalid number of parameters: must be 3");
	}
	if (typeof x !== "number" || typeof y !== "number" || typeof z !== "number") {
		throw new TypeError("Mat4Scale: Invalid parameter: must be 3 numbers");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[0][0] = x;
	this.mat[1][1] = y;
	this.mat[2][2] = z;
};
Mat4Scale.prototype = Object.create(Mat4.prototype);
Mat4Scale.prototype.constructor = Mat4Scale;
Mat4Scale.prototype.parent = Mat4;

/**
 * Creates a transformation matrix 4×4 to translate in 3D
 * @augments {Mat4}
 * @param {number} x translation along X-axis
 * @param {number} y translation along Y-axis
 * @param {number} z translation along Z-axis
 * @throws {TypeError} If any of the parameters is not defined or is not a number
 */
const Mat4Transl = function(x, y, z) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4Transl: Invalid number of parameters: must be 3");
	}
	if (typeof x !== "number" || typeof y !== "number" || typeof z !== "number") {
		throw new TypeError("Mat4Transl: Invalid parameter: must be 3 numbers");
	}

	this.mat = new Mat4Identity().mat;
	this.mat[3][0] = x;
	this.mat[3][1] = y;
	this.mat[3][2] = z;
};
Mat4Transl.prototype = Object.create(Mat4.prototype);
Mat4Transl.prototype.constructor = Mat4Transl;
Mat4Transl.prototype.parent = Mat4;

/**
 * Creates a 4×4 transition matrix from the current frame (coordinate system) to the observer (camera)
 * @augments {Mat4}
 * @param {number} e eye, position of the observer frame origin
 * @param {number} v view vector, the direction of the view of the observer
 * @param {number} u up vector
 * @throws {TypeError} If any of the parameters is not defined or is not a Vec3D
 */
const Mat4ViewRH = function(e, v, u) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4ViewRH: Invalid number of parameters: must be 3");
	}
	if (!(e instanceof Vec3D) || !(v instanceof Vec3D) || !(u instanceof Vec3D)) {
		throw new TypeError("Mat4ViewRH: Invalid parameter: must be 3 Vec3D");
	}

	this.mat = new Mat4Identity().mat;
	const z = v.mul(-1.0).normalized();
	const x = u.cross(z).normalized();
	const y = z.cross(x);

	this.mat[0][0] = x.x;
	this.mat[1][0] = x.y;
	this.mat[2][0] = x.z;
	this.mat[3][0] = -e.dot(x);
	this.mat[0][1] = y.x;
	this.mat[1][1] = y.y;
	this.mat[2][1] = y.z;
	this.mat[3][1] = -e.dot(y);
	this.mat[0][2] = z.x;
	this.mat[1][2] = z.y;
	this.mat[2][2] = z.z;
	this.mat[3][2] = -e.dot(z);
};
Mat4ViewRH.prototype = Object.create(Mat4.prototype);
Mat4ViewRH.prototype.constructor = Mat4ViewRH;
Mat4ViewRH.prototype.parent = Mat4;

/**
 * Returns the result of element-wise summation with another 4×4 matrix
 * @param {Mat3} m 4×4 matrix
 * @return {Mat4}  new Mat4 instance
 * @throws {TypeError} If m is not a Mat4
 */
Mat4.prototype.add = function(m) {
	if (m instanceof Mat4) {
		const temp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				temp.mat[i][j] = this.mat[i][j] + m.mat[i][j];
			}
		}
		return temp;
	} else {
		throw new TypeError("Mat4, add: Invalid parameter: must be a Mat4");
	}
};

/**
 * Returns the result of element-wise multiplication by the given scalar value or the result of matrix multiplication
 * @param  {number,Mat4} m a number or a Mat4
 * @return {Mat4}          new Mat4 instance
 * @throws {TypeError} If m is not a number or a Mat4
 */
Mat4.prototype.mul = function(m) {
	if (typeof m === "number") {
		const temp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				temp.mat[i][j] = this.mat[i][j] * m;
			}
		}
		return temp;
	} else if (m instanceof Mat4) {
		const temp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				let sum = 0;
				for (let k = 0; k < 4; k++) {
					sum += this.mat[i][k] * m.mat[k][j];
				}
				temp.mat[i][j] = sum;
			}
		}
		return temp;
	} else {
		throw new TypeError("Mat4, mul: Invalid parameter: must be a number or a Mat4");
	}
};

/**
 * Returns the transposition of this 4×4 matrix
 * @return {Mat4} new Mat4 instance
 */
Mat4.prototype.transpose = function() {
	const temp = new Mat4();
	for (let i = 0; i < 4; i++) {
		for (let j = 0; j < 4; j++) {
			temp.mat[i][j] = this.mat[j][i];
		}
	}
	return temp;
};

/**
 * Convert this 4×4 matrix to 3×3 by ignoring last row and last column
 * @return {Mat3} new Mat3 instance
 */
Mat4.prototype.toMat3 = function() {
	const temp = new Mat3();
	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			temp.mat[i][j] = this.mat[i][j];
		}
	}
	return temp;
};

/**
 * Print values to the console
 * @return {Mat4} reference to this instance for calls chaining
 */
Mat4.prototype.c = function() {
	window.console.log("%cMat4", 'font-style: italic; font-weight: bold');
	for (let i = 0; i < 4; i++) {
		let x = "";
		for (let j = 0; j < 4; j++) {
			x += this.mat[i][j] + ", ";
		}
		window.console.log(x.substring(0, x.length - 2));
	}
	return this;
};


/**
 * Virtual camera, controls view transformation via observer position, azimuth and zenith (in radians)
 * @constructor
 */
const Camera = function() {
	this.azimuth = 0.0;
	this.zenith = 0.0;
	this.radius = 1.0;
	this.xy = true;
	this.pos = new Vec3D(0.0, 0.0, 0.0);
	this.firstPerson = true; // true -> first person view, false -> third person view

	this.eye = new Vec3D();
	this.eyeVector = new Vec3D();
	this.view = new Mat4();

	this.computeMatrix();
};

/**
 * Recalculation of eye, eyeVector and view
 */
Camera.prototype.computeMatrix = function() {
	this.eyeVector = new Vec3D(
		Math.sin(-this.azimuth) * Math.cos(this.zenith),
		Math.cos(-this.azimuth) * Math.cos(this.zenith),
		Math.sin(this.zenith)
	);
	if (this.firstPerson) {
		this.eye = new Vec3D(this.pos);
		this.view = new Mat4ViewRH(
			this.pos,
			this.eyeVector.mul(this.radius),
			new Vec3D(
				Math.sin(-this.azimuth) * Math.cos(this.zenith + Math.PI / 2),
				Math.cos(-this.azimuth) * Math.cos(this.zenith + Math.PI / 2),
				Math.sin(this.zenith + Math.PI / 2)
			)
		);
	} else {
		this.eye = this.pos.add(this.eyeVector.mul(-1 * this.radius));
		this.view = new Mat4ViewRH(
			this.eye,
			this.eyeVector.mul(this.radius),
			new Vec3D(
				Math.sin(-this.azimuth) * Math.cos(this.zenith + Math.PI / 2),
				Math.cos(-this.azimuth) * Math.cos(this.zenith + Math.PI / 2),
				Math.sin(this.zenith + Math.PI / 2)
			)
		);
	}
};

/**
 * Adds value to current azimuth value
 * @param {number} value
 * @throws {TypeError} If value is not a number
 */
Camera.prototype.addAzimuth = function(value) {
	if (typeof value === "number") {
		this.azimuth += value;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, addAzimuth: Invalid parameter: must be a number");
	}
};

/**
 * Set azimuth to a new value
 * @param {number} value
 * @throws {TypeError} If value is not a number
 */
Camera.prototype.setAzimuth = function(value) {
	if (typeof value === "number") {
		this.azimuth = value;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setAzimuth: Invalid parameter: must be a number");
	}
};

/**
 * Adds value to current zenith value
 * @param {number} value
 * @throws {TypeError} If value is not a number
 */
Camera.prototype.addZenith = function(value) {
	if (typeof value === "number") {
		if (Math.abs(this.zenith + value) <= Math.PI / 2) {
			this.zenith += value;
			this.computeMatrix();
		}
	} else {
		throw new TypeError("Camera, addZenith: Invalid parameter: must be a number");
	}
};

/**
 * Set zenith to a new value
 * @param {number} value
 * @throws {TypeError} If value is not a number
 */
Camera.prototype.setZenith = function(value) {
	if (typeof value === "number") {
		this.zenith = value;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setZenith: Invalid parameter: must be a number");
	}
};

/**
 * Adds value to current radius value
 * @param {number} dist
 * @throws {TypeError} If dist is not a number
 */
Camera.prototype.addRadius = function(dist) {
	if (typeof dist === "number") {
		if (this.radius + dist < 0.1) {
			this.radius = 0.1;
		} else {
			this.radius += dist;
			this.computeMatrix();
		}
	} else {
		throw new TypeError("Camera, addRadius: Invalid parameter: must be a number");
	}
};

/**
 * Set radius to a new value
 * @param {number} dist
 * @throws {TypeError} If dist is not a number
 */
Camera.prototype.setRadius = function(dist) {
	if (typeof dist === "number") {
		this.radius = dist;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setRadius: Invalid parameter: must be a number");
	}
};

/**
 * Multiply zenith by given value
 * @param  {number} scale
 * @throws {TypeError} If scale is not a number
 */
Camera.prototype.mulRadius = function(scale) {
	if (typeof scale === "number") {
		if (this.radius * scale < 0.1) {
			this.radius = 0.1;
		} else {
			this.radius *= scale;
			this.computeMatrix();
		}
	} else {
		throw new TypeError("Camera, mulRadius: Invalid parameter: must be a number");
	}
};

/**
 * Move in the direction of the view vector by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.forward = function(distance) {
	if (typeof distance === "number") {
		if (!this.xy) {
			this.pos = this.pos.add(new Vec3D(
				(Math.cos(this.azimuth - Math.PI / 2) * Math.cos(this.zenith + Math.PI)),
				(Math.sin(this.azimuth - Math.PI / 2) * Math.cos(this.zenith + Math.PI)),
				Math.sin(this.zenith)).mul(distance)
			);
		} else {
			this.pos = this.pos.add(
				new Vec3D(
					Math.cos(this.azimuth - Math.PI / 2),
					Math.sin(this.azimuth - Math.PI / 2),
					0.0
				).mul(-distance)
			);
		}
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, forward: Invalid parameter: must be a number");
	}
};

/**
 * Move in the opposite direction of the view vector by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.backward = function(distance) {
	if (typeof distance === "number") {
		this.forward(-distance);
	} else {
		throw new TypeError("Camera, backward: Invalid parameter: must be a number");
	}
};

/**
 * Move to the right from the observer's perspective by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.right = function(distance) {
	if (typeof distance === "number") {
		this.pos = this.pos.add(
			new Vec3D(Math.cos(this.azimuth), Math.sin(this.azimuth), 0).mul(distance)
		);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, right: Invalid parameter: must be a number");
	}
};

/**
 * Move to the left from the observer's perspective by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.left = function(distance) {
	if (typeof distance === "number") {
		this.right(-distance);
	} else {
		throw new TypeError("Camera, left: Invalid parameter: must be a number");
	}
};

/**
 * Move in the negative direction on the Z-axis by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.down = function(distance) {
	if (typeof distance === "number") {
		this.pos.z -= distance;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, down: Invalid parameter: must be a number");
	}
};

/**
 * Move in the positive direction on the Z-axis by the given distance
 * @param  {number} distance
 * @throws {TypeError} If distance is not a number
 */
Camera.prototype.up = function(distance) {
	if (typeof distance === "number") {
		this.pos.z += distance;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, up: Invalid parameter: must be a number");
	}
};

/**
 * Move in the direction given by the input vector Vec3D
 * @param  {Vec3D} direction
 * @throws {TypeError} If direction is not a Vec3D
 */
Camera.prototype.move = function(direction) {
	if (direction instanceof Vec3D) {
		this.pos = this.pos.add(direction);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, move: Invalid parameter: must be a Vec3D");
	}
};

/**
 * Set a new position to the camera
 * @param  {Vec3D} position
 * @throws {TypeError} If position is not a Vec3D
 */
Camera.prototype.setPosition = function(position) {
	if (position instanceof Vec3D) {
		this.pos = new Vec3D(position);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setPosition: Invalid parameter: must be a Vec3D");
	}
};

/**
 * Set first/third person parameter
 * @param  {boolean} value true -> first person; false-> third person
 * @throws {TypeError} If value is not a boolean
 */
Camera.prototype.setFirstPerson = function(value) {
	if (typeof value === "boolean") {
		this.firstPerson = value;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setFirstPerson: Invalid parameter: must be a boolean");
	}
};

/**
 * Object for working with colors
 * @param {number,Col,Point3D} ar red channel (integer int <0; 255>, double in <0; 1>) or another Col or a Point3D
 * @param {number} ag             green channel (integer int <0; 255>, double in <0; 1>)
 * @param {number} ab             blue channel (integer int <0; 255>, double in <0; 1>)
 * @param {number} aa             alpha channel (integer int <0; 255>, double in <0; 1>)
 * @constructor
 * @throws {TypeError} If the parameters are not passing any of the conditions
 */
const Col = function (ar, ag, ab, aa) {
	// copy another Col
	if (ar instanceof Col) {
		this.r = ar.r;
		this.g = ar.g;
		this.b = ar.b;
		this.a = ar.a;
	// make the coordinates of a Point3D color channels
	} else if (ar instanceof Point3D) {
		this.r = ar.x;
		this.g = ar.y;
		this.b = ar.z;
		this.a = ar.w;
	// if there are 3 or 4 parameters and they are integers then it is assumed interval <0; 255>
	// !!! warning: for example r=1, g=0, b=1 passes the condition and is treated as in <0; 255> interval
	} else if ((arguments.length === 3 || arguments.length === 4) && this.isInt(ar) && this.isInt(ag) && this.isInt(ab)) {
		this.r = ar / 255.0;
		this.g = ag / 255.0;
		this.b = ab / 255.0;
		this.a = (this.isInt(aa)) ? aa / 255.0 : 1.0;
	// if there are 3 or 4 parameters and they are not integers then it is assumed interval <0; 1>
	// assumed 1.0 if the parameter is not set
	} else if (arguments.length === 3 || arguments.length === 4) {
		this.r = (typeof ar === "number") ? ar : 1.0;
		this.g = (typeof ag === "number") ? ag : 1.0;
		this.b = (typeof ab === "number") ? ab : 1.0;
		this.a = (typeof aa !== "undefined" && typeof aa === "number") ? aa : 1.0;
	} else {
		throw new TypeError("Col: Invalid parameter.");
	}
};

/**
 * Check if a number is an integer
 * @param  {number}  n number to check
 * @return {Boolean}   true if number is an integer, false otherwise
 */
Col.prototype.isInt = function(n) {
	return (n % 1 === 0);
};

/**
 * Add another Col to this Col, ignore alpha channel
 * @param  {Col} c Col to add
 * @return {Col}   new Col instance
 * @throws {TypeError} If c is not a Col
 */
Col.prototype.addna = function(c) {
	if (c instanceof Col) {
		return new Col(this.r + c.r, this.g + c.g, this.b + c.b);
	} else {
		throw new TypeError("Col, addna: Invalid parameter: must be a Col");
	}
};

/**
 * Multiply this Col by a number, ignore alpha channel
 * @param  {number} x number
 * @return {Col}      new Col instance
 * @throws {TypeError} If x is not a number
 */
Col.prototype.mulna = function(x) {
	if (typeof x === "number") {
		return new Col(this.r * x, this.g * x, this.b * x);
	} else {
		throw new TypeError("Col, mulna: Invalid parameter: must be a number");
	}
};

/**
 * Add another Col to this Col, including alpha channel
 * @param {Col} c Col to add
 * @return {Col}  new Col instance
 * @throws {TypeError} If c is not a Col
 */
Col.prototype.add = function(c) {
	if (c instanceof Col) {
		return new Col(this.r + c.r, this.g + c.g, this.b + c.b, this.a + c.a);
	} else {
		throw new TypeError("Col, add: Invalid parameter: must be a Col");
	}
};

/**
 * Multiply this Col by a number or another Col, including alpha channel
 * @param  {number,Col} c a number or another Col
 * @return {Col}          new Col instance
 * @throws {TypeError} If c is not a Col or a number
 */
Col.prototype.mul = function(c) {
	if (c instanceof Col) {
		return new Col(this.r * c.r, this.g * c.g, this.b * c.b, this.a * c.a);
	} else if (typeof c === "number") {
		return new Col(this.r * c, this.g * c, this.b * c, this.a * c);
	} else {
		throw new TypeError("Col, mul: Invalid parameter: must be a Col or a number");
	}
};

/**
 * Returns a new Col with gamma-correction applied to RGB channels
 * @param  {number} gamma gamma value
 * @return {Col}          new Col instance
 * @throws {TypeError} If gamma is not a number
 */
Col.prototype.gamma = function(gamma) {
	if (typeof gamma === "number") {
		return new Col(
			Math.pow(this.r, gamma),
			Math.pow(this.g, gamma),
			Math.pow(this.b, gamma),
			this.a
		);
	} else {
		throw new TypeError("Col, gamma: Invalid parameter: must be a number");
	}
};

/**
 * Returns a new Col with clamped channels to <0;1> (except for alpha)
 * @return {Col} new Col instance
 */
Col.prototype.saturate = function() {
	return new Col(
		Math.max(0, Math.min(this.r, 1)),
		Math.max(0, Math.min(this.g, 1)),
		Math.max(0, Math.min(this.b, 1)),
		this.a
	);
};

/**
 * Returns the RGB channels scaled to <0;255> and packed to individual bytes of a 32-bit integer value
 * with blue in the least significant byte and zero in the most significant byte.
 * @return {number} number with RGB in lower 24 bits of a 32-bit integer
 */
Col.prototype.getRGB = function() {
	return ((this.r * 255.0) << 16) | ((this.g * 255.0) << 8) | (this.b * 255.0);
};

/**
 * Returns the ARGB channels scaled to <0;255> and packed to individual bytes of a 32-bit integer value
 * with blue in the least significant byte and alpha in the most significant byte.
 * @return {number} number with ARGB in a 32-bit integer
 */
Col.prototype.getARGB = function() {
	return ((this.a * 255.0) << 24) | ((this.r * 255.0) << 16) | ((this.g * 255.0) << 8) | (this.b * 255.0);
};

/**
 * Print values to the console
 * @return {Col} reference to this instance for calls chaining
 */
Col.prototype.c = function() {
	window.console.log(this);
	return this;
};


/**
 * Object for working with cubics - Bézier, Ferguson and coons
 * @param {number} type 1 -> Ferguson, 2 -> coons, other -> Bézier
 * @constructor
 */
const Kubika = function(type) {
	// base matrix
	this.bm = new Mat4();
	// control points matrix
	this.rb;

	switch (typeof type === "number" ? type : 0) {
		case 1: // Ferguson
			this.type = 1;

			this.bm.mat[0][0] = 2;
			this.bm.mat[0][1] = -2;
			this.bm.mat[0][2] = 1;
			this.bm.mat[0][3] = 1;

			this.bm.mat[1][0] = -3;
			this.bm.mat[1][1] = 3;
			this.bm.mat[1][2] = -2;
			this.bm.mat[1][3] = -1;

			this.bm.mat[2][0] = 0;
			this.bm.mat[2][1] = 0;
			this.bm.mat[2][2] = 1;
			this.bm.mat[2][3] = 0;

			this.bm.mat[3][0] = 1;
			this.bm.mat[3][1] = 0;
			this.bm.mat[3][2] = 0;
			this.bm.mat[3][3] = 0;
			break;

		case 2: // Coons
			this.type = 2;

			this.bm.mat[0][0] = -1;
			this.bm.mat[0][1] = 3;
			this.bm.mat[0][2] = -3;
			this.bm.mat[0][3] = 1;

			this.bm.mat[1][0] = 3;
			this.bm.mat[1][1] = -6;
			this.bm.mat[1][2] = 3;
			this.bm.mat[1][3] = 0;

			this.bm.mat[2][0] = -3;
			this.bm.mat[2][1] = 0;
			this.bm.mat[2][2] = 3;
			this.bm.mat[2][3] = 0;

			this.bm.mat[3][0] = 1;
			this.bm.mat[3][1] = 4;
			this.bm.mat[3][2] = 1;
			this.bm.mat[3][3] = 0;

			this.bm = this.bm.mul(1 / 6);
			break;

		case 0: // Bézier
		default:
			this.type = 0;

			this.bm.mat[0][0] = -1;
			this.bm.mat[0][1] = 3;
			this.bm.mat[0][2] = -3;
			this.bm.mat[0][3] = 1;

			this.bm.mat[1][0] = 3;
			this.bm.mat[1][1] = -6;
			this.bm.mat[1][2] = 3;
			this.bm.mat[1][3] = 0;

			this.bm.mat[2][0] = -3;
			this.bm.mat[2][1] = 3;
			this.bm.mat[2][2] = 0;
			this.bm.mat[2][3] = 0;

			this.bm.mat[3][0] = 1;
			this.bm.mat[3][1] = 0;
			this.bm.mat[3][2] = 0;
			this.bm.mat[3][3] = 0;
	}
};

/**
 * Initialize with 4 control points
 * @param  {Point3D} b1 control point
 * @param  {Point3D} b2 control point
 * @param  {Point3D} b3 control point
 * @param  {Point3D} b4 control point
 * @throws {TypeError} If any of the parameters is not set or is not a Point3D
 */
Kubika.prototype.init = function(b1, b2, b3, b4) {
	if (arguments.length !== 4) {
		throw new TypeError("Kubika, init: Invalid number of parameters: must be 4");
	} else if (b1 instanceof Point3D && b2 instanceof Point3D && b3 instanceof Point3D && b4 instanceof Point3D) {
		if (this.type === 1) {
			this.rb = new Mat4(b1, b4, b2.sub(b1), b4.sub(b3));
		} else {
			this.rb = new Mat4(b1, b2, b3, b4);
		}
		this.rb = this.bm.mul(this.rb);
	} else {
		throw new TypeError("Kubika, init: Invalid parameter: must be 4 Point3D");
	}
};

/**
 * Compute the coordinates of a point on the cubic curve corresponding to the parameter from <0; 1>
 * @param  {number} t parameter from interval <0; 1>, clamped is outside of this range
 * @return {Point3D}  new Point3D on the curve
 * @throws {TypeError} If t is not a number
 */
Kubika.prototype.compute = function(t) {
	if (typeof t === "number") {
		if (t > 1) t = 1;
		if (t < 0) t = 0;

		const res1 = new Point3D(t * t * t, t * t, t, 1);
		const res2 = res1.mul(this.rb);
		res2.w = 1;
		return res2;
	} else {
		throw new TypeError("Kubika, compute: Invalid parameter: must be a number");
	}
};

/**
 * Object for working with bicubics - Bézier, Ferguson and Coons
 * @param {number} type 1 -> Ferguson, 2 -> coons, other -> Beziér
 * @constructor
 */
const Bikubika = function(type) {
	// Point3D
	this.u1; this.u2; this.u3; this.u4;

	this.k1 = new Kubika(type);
	this.k2 = new Kubika(type);
	this.k3 = new Kubika(type);
	this.k4 = new Kubika(type);
	this.k5 = new Kubika(type);
};

/**
 * Initialization with 4×4 control points
 * @param  {Point3D} b11 control point
 * @param  {Point3D} b12 control point
 * @param  {Point3D} b13 control point
 * @param  {Point3D} b14 control point
 * @param  {Point3D} b21 control point
 * @param  {Point3D} b22 control point
 * @param  {Point3D} b23 control point
 * @param  {Point3D} b24 control point
 * @param  {Point3D} b31 control point
 * @param  {Point3D} b32 control point
 * @param  {Point3D} b33 control point
 * @param  {Point3D} b34 control point
 * @param  {Point3D} b41 control point
 * @param  {Point3D} b42 control point
 * @param  {Point3D} b43 control point
 * @param  {Point3D} b44 control point
 */
Bikubika.prototype.init = function(
			b11, b12, b13, b14, b21, b22, b23, b24,
			b31, b32, b33, b34, b41, b42, b43, b44) {
	this.k1.init(b11, b12, b13, b14);
	this.k2.init(b21, b22, b23, b24);
	this.k3.init(b31, b32, b33, b34);
	this.k4.init(b41, b42, b43, b44);
};

/**
 * Compute the coordinates of a point on the bicubic surface corresponding to the u,v parameters from <0; 1>
 * @param  {number} u parameter from interval <0; 1>, clamped is outside of this range
 * @param  {number} v parameter from interval <0; 1>, clamped is outside of this range
 * @return {Point3D}  new Point3D on the surface
 * @throws {TypeError} If u and v are not numbers
 */
Bikubika.prototype.compute = function(u, v) {
	if (typeof u === "number" && typeof v === "number") {
		if (u > 1) u = 1;
		if (u < 0) u = 0;
		if (v > 1) v = 1;
		if (v < 0) v = 0;

		this.u1 = this.k1.compute(u);
		this.u2 = this.k2.compute(u);
		this.u3 = this.k3.compute(u);
		this.u4 = this.k4.compute(u);
		this.k5.init(this.u1, this.u2, this.u3, this.u4);
		return this.k5.compute(v);
	} else {
		throw new TypeError("Bikubika, compute: Invalid parameter: must be 2 numbers");
	}
};
