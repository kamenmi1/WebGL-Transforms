"use strict";

/**
 * @description Tento soubor obsahuje soubor objektů a funkcí umožňující
 *  použivat vektorové a maticové výpočty v JavaScriptu. Navíc obsahuje
 *  funcionalitu pro práci s kvaterniony, kubikami a bikubikami
 *  a speciální objekt pro práci s kamerou.
 *
 * 	Založeno na balíčku javovských tříd "transforms3d" použivaného na FIM UHK při výuce předmětu PGRF.
 * 	Největší změny oproti vzoru:
 * 		- ZeroArray: v JavaScriptu je odlišná práce s poli, nutný pomocný objekt
 * 		- Mat4: toMat3 metoda
 * 		- Camera: změny pro uzpůsobení WebGL.
 * 		- Camera: přidána vlastnost "xy" pro možnost uzamknutí pohybu v rámci roviny os x a y
 * 		- Col: vynechány některé konstruktory
 * 		- "c" funkce pro výpis do konzole
 * 
 * @version 1.1
 * (June 2017)
 */

/**
 * Vytvoření matice o zvoleném rozměru obsahující nuly
 * @param {number} x rozměr matice (např. 3 pro 3×3 matici), pokud není zadáno, tak se předpokládá matice 4×4
 * @throws {TypeError} If x je zadáno a není číslo
 * @constructor
 */
var ZeroArray = function(x) {
	let length = (typeof x !== "undefined") ? x : 4;
	if (typeof length === "number") {
		let mat = [];
		for (let i = 0; i < length; i++) {
			mat[i] = [];
			for (let j = 0; j < length; j++) {
				mat[i][j] = 0.0;
			}
		}
		return mat;
	} else {
		throw new TypeError("ZeroArray: Neplatný parametr: musí být číslo");
	}
};

/**
 * Objekt pro práci s 1D vektory
 * @param {Vec1D,number} x volitelné, pokud není zadáno, tak 0
 * @constructor
 */
var Vec1D = function(x) {
	if (x instanceof Vec1D) {
		this.x = x.x;
	} else {
		if (typeof x !== "undefined") {
			if (typeof x === "number") {
				this.x = x;
			} else {
				throw new TypeError("Vec1D: Neplatný parametr x: musí být Vec1D nebo number")
			}
		} else {
			this.x = 0.0;
		}
	}
};

/**
 * Přičtení vektoru
 * @param {Vec1D} v vektor pro přičtení
 * @return {Vec1D}  nová instance Vec1D
 * @throws {TypeError} If v není Vec1D
 */
Vec1D.prototype.add = function(v) {
	if (v instanceof Vec1D) {
		return new Vec1D(this.x + v.x);
	} else {
		throw new TypeError("Vec1D, add: Neplatný parametr: musí být Vec1D");
	}
};

/**
 * Násobení skalárem
 * @param  {number} m skalár
 * @return {Vec1D}    nová instance Vec1D
 * @throws {TypeError} If m není číslo
 */
Vec1D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec1D(this.x * m);
	} else {
		throw new TypeError("Vec1D, mul: Neplatný parametr: musí být číslo");
	}
};

/**
 * Výpis hodnoty do konzole
 * @return {Vec1D} reference na volanou instanci pro možné další řetězení
 */
Vec1D.prototype.c = function() {
	console.log(this);
	return this;
};

/**
 * Objekt pro práci s 2D vektory
 * @param {number,Vec1D} x volitelné, pokud není, tak 0
 * @param {number} y       volitelné, pokud není, tak 0
 */
var Vec2D = function(x, y) {
	if (x instanceof Vec2D) {
		this.x = x.x;
		this.y = x.y;
	} else {
		// x
		if (typeof x !== "undefined") {
			if (typeof x === "number") {
				this.x = x;
			} else {
				throw new TypeError("Vec2D: Neplatný parametr x: musí být Vec2D nebo number")
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof y !== "undefined") {
			if (typeof y === "number") {
				this.y = y;
			} else {
				throw new TypeError("Vec2D: Neplatný parametr y: musí být Vec2D nebo number")
			}
		} else {
			this.y = 0.0;
		}
	}
};

/**
 * Přičtení vektoru
 * @param {Vec2D} v vektor (x,y)
 * @return {Vec2D}  nová instance Vec2D
 * @throws {TypeError} If v není Vec2D
 */
Vec2D.prototype.add = function(v) {
	if (v instanceof Vec2D) {
		return new Vec2D(this.x + v.x, this.y + v.y);
	} else {
		throw new TypeError("Vec2D, add: Neplatný parametr: musí být Vec2D");
	}
};

/**
 * Násobení skalárem, násobení vektorem
 * @param  {number,Vec2D} m skalár, vektor (x, y)
 * @return {Vec2D}          nová instance Vec2D
 * @throws {TypeError} If m nená číslo ani Vec2D
 */
Vec2D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec2D(this.x * m, this.y * m);
	} else if (m instanceof Vec2D) {
		return new Vec2D(this.x * m.x, this.y * m.y);
	} else {
		throw new TypeError("Vec2D, mul: Neplatný parametr: musí být číslo nebo Vec2D");
	}
};

/**
 * Skalární součin vektoru
 * @param  {Vec2D} v vektor (x, y)
 * @return {Vec2D}   nová instance Vec2D
 * @throws {TypeError} If v není Vec2D
 */
Vec2D.prototype.dot = function(v) {
	if (v instanceof Vec2D) {
		return this.x * v.x + this.y * v.y;
	} else {
		throw new TypeError("Vec2D, dot: Neplatný parametr: musí být Vec2D");
	}
};

/**
 * Velikost vektoru
 * @return {number} velikost
 */
Vec2D.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y);
};

/**
 * Normalizace vektoru
 * @return {Vec2D} nová instance Vec2D
 */
Vec2D.prototype.normalized = function() {
	let len = this.length();
	if (len === 0.0) return new Vec2D(0, 0);
	return new Vec2D(this.x / len, this.y / len);
};

/**
 * Výpis hodnoty do konzole
 * @return {Vec2D} reference na volanou instanci pro možné další řetězení
 */
Vec2D.prototype.c = function() {
	console.log(this);
	return this;
};


/**
 * Objekt pro práci s 3D vektory
 * @param {number,Vec3D,Point3D} x volitelný agrument, pokud není zadán, tak 0
 * @param {number} y volitelný agrument, pokud není zadán, tak 0
 * @param {number} z volitelný agrument, pokud není zadán, tak 0
 * @constructor
 */
var Vec3D = function(x, y, z) {
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
				throw new TypeError("Vec3D: Neplatný parametr x: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof y !== "undefined") {
			if (typeof y === "number") {
				this.y = y;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr y: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.y = 0.0;
		}
		// z
		if (typeof z !== "undefined") {
			if (typeof z === "number") {
				this.z = z;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr z: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.z = 0.0;
		}
	}
};

/**
 * Přičtení vektoru
 * @param {Vec3D} v vektor (x, y, z)
 * @return {Vec3D} nová instance Vec3D
 * @throws {TypeError} If v není Vec3D
 */
Vec3D.prototype.add = function(v) {
	if (v instanceof Vec3D) {
		return new Vec3D(this.x + v.x, this.y + v.y, this.z + v.z);
	} else {
		throw new TypeError("Vec3D, add: Neplatný parametr: musí být Vec3D");
	}
};

/**
 * Násobení skalárem, maticí zprava, kvaterninonem, vektoem po složkách
 * @param  {number,Mat3,Vec3D,Quat} m skalár, matice3×3, kvaternion, vektor (x, y, z)
 * @return {Vec3D}   nová instance Vec3D
 * @throws {TypeError} If m nemá povolený typ
 */
Vec3D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Vec3D(this.x * m, this.y * m, this.z * m);
	} else if (m instanceof Mat3) {
		let res = new Vec3D();
		res.x = m.mat[0][0] * this.x + m.mat[1][0] * this.y + m.mat[2][0] * this.z;
		res.y = m.mat[0][1] * this.x + m.mat[1][1] * this.y + m.mat[2][1] * this.z;
		res.z = m.mat[0][2] * this.x + m.mat[1][2] * this.y + m.mat[2][2] * this.z;
		return res;
	} else if (m instanceof Vec3D) {
		return new Vec3D(this.x * m.x, this.y * m.y, this.z * m.z);
	} else if (m instanceof Quat) {
		let p = new Quat(0, this.x, this.y, this.z);
		p = m.mulR(p).mulR(m.inv());
		return new Vec3D(p.i, p.j, p.k);
	} else {
		throw new TypeError("Vec3D, mul: Neplatný parametr: musí být číslo, Mat3, Vec3D nebo Quat");
	}
};

/**
 * Skalární součin vektoru
 * @param  {Vec3D} rhs vektor (x, y, z)
 * @return {number}    součin
 * @throws {TypeError} If rhs není Vec3D
 */
Vec3D.prototype.dot = function(rhs) {
	if (rhs instanceof Vec3D) {
		return this.x * rhs.x + this.y * rhs.y + this.z * rhs.z;
	} else {
		throw new TypeError("Vec3D, dot: Neplatný parametr: musí být Vec3D");
	}
};

/**
 * Vektorový součin vektorů
 * @param  {Vec3D} v vektor (x, y, z)
 * @return {Vec3D}   nová instance Vec3D
 * @throws {TypeError} If v není Vec3D
 */
Vec3D.prototype.cross = function(v) {
	if (v instanceof Vec3D) {
		return new Vec3D(this.y * v.z - this.z * v.y, this.z * v.x - this.x * v.z, this.x * v.y - this.y * v.x);
	} else {
		throw new TypeError("Vec3D, cross: Neplatný parametr: musí být Vec3D");
	}
};

/**
 * Veliskot vektoru
 * @return {number} velikost
 */
Vec3D.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
};

/**
 * Normalizace vektoru
 * @return {Vec3D} nová instance Vec3D
 */
Vec3D.prototype.normalized = function() {
	let len = this.length();
	if (len === 0) return new Vec3D(0, 0, 0);
	return new Vec3D(this.x / len, this.y / len, this.z / len);
};

/**
 * Výpis hodnoty do konzole
 * @return {Vec3D} reference na volanou instanci pro možné další řetězení
 */
Vec3D.prototype.c = function() {
	console.log(this);
	return this;
};


/**
 * Objekt pro práci s body ve 3D (homogenní souřadnice)
 * @param {number,Vec3D,Point3D} ax volitelný argument, pokud není zadán, tak 0
 * @param {number} ay volitelný argument, pokud není zadán, tak 0
 * @param {number} az volitelný argument, pokud není zadán, tak 0
 * @param {number} aw volitelný argument, pokud není zadán, tak 0
 * @constructor
 */
var Point3D = function(ax, ay, az, aw) {
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
				throw new TypeError("Vec3D: Neplatný parametr ax: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.x = 0.0;
		}
		// y
		if (typeof ay !== "undefined") {
			if (typeof ay === "number") {
				this.y = ay;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr y: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.y = 0.0;
		}
		// z
		if (typeof az !== "undefined") {
			if (typeof az === "number") {
				this.z = az;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr z: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.z = 0.0;
		}
		// w
		if (typeof aw !== "undefined") {
			if (typeof aw === "number") {
				this.w = aw;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr w: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.w = 1.0;
		}
	}
};

/**
 * Násobení skalárem, násobení maticí zprava, transformace 3D bodu kvaternionem
 * @param  {number,Mat4,Quat} m skalár, matice 4×4, kvaternion
 * @return {Point3D}   nová instance Point3D
 * @throws {TypeError} If m není číslo, Mat4 ani Quat
 */
Point3D.prototype.mul = function(m) {
	if (typeof m === "number") {
		return new Point3D(this.x * m, this.y * m, this.z * m, this.w * m);
	} else if (m instanceof Mat4) {
		let res = new Point3D();
		res.x = m.mat[0][0] * this.x + m.mat[1][0] * this.y + m.mat[2][0] * this.z + m.mat[3][0] * this.w;
		res.y = m.mat[0][1] * this.x + m.mat[1][1] * this.y + m.mat[2][1] * this.z + m.mat[3][1] * this.w;
		res.z = m.mat[0][2] * this.x + m.mat[1][2] * this.y + m.mat[2][2] * this.z + m.mat[3][2] * this.w;
		res.w = m.mat[0][3] * this.x + m.mat[1][3] * this.y + m.mat[2][3] * this.z + m.mat[3][3] * this.w;
		return res;
	} else if (m instanceof Quat) {
		return new Point3D(this.dehomog().mul(m));
	} else {
		throw new TypeError("Point3D, mul: Neplatný parametr: musí být číslo, Mat4 nebo Quat");
	}
};

/**
 * Přičtení vektoru
 * @param {Point3D} p vektor (x, y, z, w)
 * @return {Point3D} nová instance Point3D
 * @throws {TypeError} If p není Point3D
 */
Point3D.prototype.add = function(p) {
	if (p instanceof Point3D) {
		return new Point3D(this.x + p.x, this.y + p.y, this.z + p.z, this.w + p.w);
	} else {
		throw new TypeError("Point3D, add: Neplatný parametr: musí být Point3D");
	}
};

/**
 * Odečtení vektoru
 * @param  {Point3D} p vektor (x, y, z, w)
 * @return {Point3D}   nová instance Point3D
 * @throws {TypeError} If p není Point3D
 */
Point3D.prototype.sub = function(p) {
	if (p instanceof Point3D) {
		return new Point3D(this.x - p.x, this.y - p.y, this.z - p.z, this.w - p.w);
	} else {
		throw new TypeError("Point3D, sub: Neplatný parametr: musí být Point3D");
	}
};

/**
 * Dehmogenizace vektoru
 * @return {Vec3D} nová instance Vec3D
 */
Point3D.prototype.dehomog = function() {
	if (this.w === 0) return new Vec3D(0, 0, 0);
	return new Vec3D(this.x / this.w, this.y / this.w, this.z / this.w);
};

/**
 * Převod na vektor (x, y, z), zanedbání w
 * @return {Vec3D} nová instance Vec3D
 */
Point3D.prototype.ignoreW = function() {
	return new Vec3D(this.x, this.y, this.z);
};

/**
 * Výpis hodnoty do konzole
 * @return {Point3D} reference na volanou instanci pro možné další řetězení
 */
Point3D.prototype.c = function() {
	console.log(this);
	return this;
};

/**
 * Objekt pro práci s kvaterniony
 * @param {number,Quat} r volitelný argument, pokud není zadán, tak 0
 * @param {number} i volitelný argument, pokud není zadán, tak 0
 * @param {number} j volitelný argument, pokud není zadán, tak 0
 * @param {number} k volitelný argument, pokud není zadán, tak 0
 * @constructor
 */
var Quat = function(r, i, j, k) {
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
				throw new TypeError("Vec3D: Neplatný parametr i: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.i = 0.0;
		}
		// j
		if (typeof j !== "undefined") {
			if (typeof j === "number") {
				this.j = j;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr j: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.j = 0.0;
		}
		// k
		if (typeof k !== "undefined") {
			if (typeof k === "number") {
				this.k = k;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr k: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.k = 0.0;
		}
		// r
		if (typeof r !== "undefined") {
			if (typeof r === "number") {
				this.r = r;
			} else {
				throw new TypeError("Vec3D: Neplatný parametr r: musí být Vec3D, Point3D nebo number")
			}
		} else {
			this.r = 0.0;
		}
	}
};

/**
 * Sčítání kvaternionů
 * @param {Quat} q kvaternion
 * @return {Quat}  nová instance Quat
 * @throws {TypeError} If q není Quat
 */
Quat.prototype.add = function(q) {
	if (q instanceof Quat) {
		return new Quat(this.r + q.r, this.i + q.i, this.j + q.j, this.k + q.k);
	} else {
		throw new TypeError("Quat, add: Neplatný parametr: musí být Quat");
	}
};

/**
 * Odčítání kvaternionů
 * @param {Quat} q kvaternion
 * @return {Quat}  nová instance Quat
 * @throws {TypeError} If q není Quat
 */
Quat.prototype.sub = function(q) {
	if (q instanceof Quat) {
		return new Quat(this.r - q.r, this.i - q.i, this.j - q.j, this.k - q.k);
	} else {
		throw new TypeError("Quat, sub: Neplatný parametr: musí být Quat");
	}
};

/**
 * Násobení kvaternionu kvaternionem zprava
 * @param {Quat} q kvaternion
 * @return {Quat}  nová instance Quat
 * @throws {TypeError} If q není Quat
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
		throw new TypeError("Quat, mulR: Neplatný parametr: musí být Quat");
	}
};

/**
 * Násobení kvaternionu kvaternionem zleva
 * @param {Quat} q kvaternion
 * @return {Quat}  nová instance Quat
 * @throws {TypeError} If q není Quat
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
		throw new TypeError("Quat, mulL: Neplatný parametr: musí být Quat");
	}
};

/**
 * Násobení kvaterionu skalárem, kvaternionem zprava
 * @param {number,Quat} q skalár, kvaternion
 * @return {Quat}  nová instance Quat
 * @throws {TypeError} If q není Quat ani číslo
 */
Quat.prototype.mul = function(q) {
	if (typeof q === "number") {
		return new Quat(q * this.r, q * this.i, q * this.j, q * this.k);
	}
	else if (q instanceof Quat) {
		return this.mulR(q);
	} else {
		throw new TypeError("Quat, mul: Neplatný parametr: musí být Quat");
	}
};

/**
 * Velikost (norma) kvaternionu
 * @return {number} velikost
 */
Quat.prototype.norma = function() {
	return Math.sqrt(this.r * this.r + this.i * this.i + this.j * this.j + this.k * this.k);
};

/**
 * Inverzní kvaternion
 * @return {Quat} nová instance Quat
 */
Quat.prototype.inv = function() {
	let norm = this.norma();
	norm = norm * norm;
	if (norm > 0) {
		return new Quat(this.r / norm, -this.i / norm, -this.j / norm, -this.k / norm);
	} else {
		return new Quat(0, 0, 0, 0);
	}
};

/**
 * Logaritmická funkce kvaternionu
 * @return {Quat} nová instance Quat
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
		let s = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
		let a = Math.atan2(s, this.r) / s;
		return new Quat(Math.log(this.norma()), a * this.i, a * this.j, a * this.k);
	}
};

/**
 * Exponenciální funkce kvaternionu
 * @return {Quat} nová instance Quat
 */
Quat.prototype.exp = function() {
	if (this.i === 0 && this.j === 0 && this.k === 0) {
		return new Quat(Math.exp(this.r), 0, 0, 0);
	} else {
		let s = Math.sqrt(this.i * this.i + this.j * this.j + this.k * this.k);
		let cos = Math.cos(s);
		s = Math.exp(this.r) * Math.sin(s) / s;
		return new Quat(Math.exp(this.r) * cos, s * this.i, s * this.j, s * this.k);
	}
};

/**
 * Opačný kvaternion
 * @return {Quat} nová instance Quat
 */
Quat.prototype.neg = function() {
	return new Quat(-this.r, -this.i, -this.j, -this.k);
};

/**
 * Skalární součin kvaternionu
 * @param  {Quat} q  kvaternion
 * @return {number}  skalár
 */
Quat.prototype.dot = function(q) {
	if (q instanceof Quat) {
		return this.i * q.i + this.j * q.j + this.k * q.k + this.r * q.r;
	} else {
		throw new TypeError("Quat, dot: Neplatný parametr: musí být Quat");
	}
};

/**
 * Normalizace kvaternionu
 * @return {Quat} nová instance Quat
 */
Quat.prototype.renorm = function() {
	let norm = this.norma();
	if (norm > 0) {
		return new Quat(this.r / norm, this.i / norm, this.j / norm, this.k / norm);
	} else {
		return new Quat(0, 0, 0, 0);
	}
};

/**
 * Vrací matici rotace na základě kvaternionu
 * @return {Mat4} nová instance Mat4
 */
Quat.prototype.toRotationMatrix = function() {
	let res = new Mat4Identity();
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
 * Pomocný statický objekt pro práci s kvaterninony
 * @static
 */
var Quat2 = {};

/**
 * Vrací kvaternion na základě matice rotace
 * @param  {Mat4} mat
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If mat není Mat4
 */
Quat2.fromRotationMatrix = function(mat) {
	if (mat instanceof Mat4) {
		let r, i, j, k;
		let diagonal = mat.mat[0][0] + mat.mat[1][1] + mat.mat[2][2];

		if (diagonal > 0.0) {
			r = (0.5 * Math.sqrt(diagonal + mat.mat[3][3]));
			i = (mat.mat[2][1] - mat.mat[1][2]) / (4 * r);
			j = (mat.mat[0][2] - mat.mat[2][0]) / (4 * r);
			k = (mat.mat[1][0] - mat.mat[0][1]) / (4 * r);
		} else {
			let indices = [1, 2, 0];
			let a = 0;

			if (mat.mat[1][1] > mat.mat[0][0]) a = 1;
			if (mat.mat[2][2] > mat.mat[a][a]) a = 2;

			let b = indices[a];
			let c = indices[b];

			diagonal = mat.mat[a][a] - mat.mat[b][b] - mat.mat[c][c] + mat.mat[3][3];
			r = (0.5 * Math.sqrt(diagonal));
			i = (mat.mat[a][b] + mat.mat[b][a]) / (4 * r);
			j = (mat.mat[a][c] + mat.mat[c][a]) / (4 * r);
			k = (mat.mat[b][c] - mat.mat[c][b]) / (4 * r);
		}
		return new Quat(r, i, j, k);
	} else {
		throw new TypeError("Quat, fromRotationMatrix: Neplatný parametr: musí být Mat4");
	}
};

/**
 * Vrací kvaternion na základě úhlu rotace kolem jednotlivých os
 * @param  {number} a úhel rotace kolem osy x
 * @param  {number} b úhel rotace kolem osy y
 * @param  {number} c úhel rotace kolem osy z
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If některý z parametrů není číslo
 */
Quat2.fromEulerAngles = function(a, b, c) {
	if (arguments.length !== 3) {
		throw new TypeError("Quat, fromEulerAngles: Neplatný počet parametrů: musí být 3");
	} else if (typeof a === "number" && typeof b === "number" && typeof c === "number") {
		let Qi = Quat2.fromEulerAngle(a, 1, 0, 0);
		let Qj = Quat2.fromEulerAngle(b, 0, 1, 0);
		let Qk = Quat2.fromEulerAngle(c, 0, 0, 1);
		return new Quat(Qk.mul(Qj).mul(Qi));
	} else {
		throw new TypeError("Quat, fromEulerAngles: Neplatný parametr: musí být číslo");
	}
};

/**
 * Vrací kvaternion na základě úhlu a osy rotace
 * @param  {number} angle úhel rotace
 * @param  {number} a     souřadnice x osy rotace
 * @param  {number} b     souřadnice y osy rotace
 * @param  {number} c     souřadnice z osy rotace
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If některý z parametrů není číslo
 */
Quat2.fromEulerAngle = function(angle, a, b, c) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, fromEulerAngle: Neplatný počet parametrů: musí být 4");
	} else if (typeof angle === "number" && typeof a === "number" && typeof b === "number" && typeof c === "number") {
		return new Quat(
			Math.cos(angle / 2),
			Math.sin(angle / 2) * a,
			Math.sin(angle / 2) * b,
			Math.sin(angle / 2) * c
		);
	} else {
		throw new TypeError("Quat, fromEulerAngle: Neplatný parametr: musí být číslo");
	}
};

/**
 * Vrací úhel a osu rotace ve formátu Point3D(uhel,osaX,osaY,osaZ)
 * @return {Point3D} nová instance Point3D
 */
Quat.prototype.toEulerAngle = function() {
	let angle = 2 * Math.acos(this.r);
	let x = this.i;
	let y = this.j;
	let z = this.k;

	let s = Math.sqrt(x * x + y * y + z * z);
	if (s < 0.0001) {
		return new Point3D(angle, 1.0, 0.0, 0.0);
	} else {
		return new Point3D(angle, (x / s), (y / s),	(z / s));
	}
};

/**
 * Lineární interpolace pomocí kvaternionu Lerp(Q1,Q2,t)=(1-t)Q1+tQ2
 * @param  {Quat} q kvaternion
 * @param  {number} t váha z intervalu <0; 1>
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If q není Quat nebo t není číslo
 */
Quat.prototype.lerp = function(q, t) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, lerp: Neplatný počet parametrů: musí být 2");
	} else if (q instanceof Quat && typeof t === "number") {
		if (t >= 1) {
			return new Quat(q);
		} else if (t <= 0) {
			return new Quat(this);
		} else {
			return new Quat((this.mul(1 - t)).add(q.mul(t)));
		}
	} else {
		throw new TypeError("Quat, lerp: Neplatný parametr: musí být Quat a číslo");
	}
};

/**
 * Sférická interpolace pomocí kvaternionu
 * @param  {Quat} q kvaternion
 * @param  {number} t váha z intervalu <0; 1>
 * @return {Point3D} nová instance Point3D
 * @throws {TypeError} If q není Quat nebo t není číslo
 */
Quat.prototype.slerp = function(q, t) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, lerp: Neplatný počet parametrů: musí být 2");
	} else if (q instanceof Quat && typeof t === "number") {
		let c = this.dot(q);
		if (c > 1.0) {
			c = 1.0;
		} else if (c < -1.0) {
			c = -1.0;
		}
		let uhel = Math.acos(c);
		if (Math.abs(uhel) < 1.0e-5) {
			return new Quat(this);
		}
		let s = 1 / Math.sin(uhel);
		if (t >= 1) {
			return new Quat(this);
		} else if (t <= 0) {
			return new Quat(q);
		} else {
			return new Quat(
				this.renorm()
					.mul(Math.sin((1 - t) * uhel) * s)
					.add(q.renorm().mul(Math.sin(t * uhel) * s))
				).renorm();
		}
	} else {
		throw new TypeError("Quat, lerp: Neplatný parametr: musí být Quat a číslo");
	}
};

/**
 * Kubická interpolace pomocí kvaternionů
 * @param  {Quat} q   kvaternion
 * @param  {Quat} q1  kvaternion
 * @param  {Quat} q2  kvaternion
 * @param  {number} t váha z intervalu <0; 1>
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If q, q1, q2 nejsou Quat nebo t není číslo
 */
Quat.prototype.squad = function(q, q1, q2, t) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, squad: Neplatný počet parametrů: musí být 4");
	} else if (q instanceof Quat && q1 instanceof Quat && q2 instanceof Quat && typeof t === "number") {
		return new Quat(this.slerp(q, t).slerp(q1.slerp(q2, t), (2 * t * (1 - t))));
	} else {
		throw new TypeError("Quat, squad: Neplatný parametr: musí být Quat a číslo");
	}
};

/**
 * 
 * @param  {Quat} q1 kvaternion
 * @param  {Quat} q2 kvaternion
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If q1 nebo q2 není Quat
 */
Quat.prototype.quadrangle = function(q1, q2) {
	if (arguments.length !== 2) {
		throw new TypeError("Quat, quadrangle: Neplatný počet parametrů: musí být 2");
	} else if (q1 instanceof Quat && q2 instanceof Quat) {
		let s1 = this.inv().mul(q1);
		let s2 = this.inv().mul(q2);
		return new Quat((s1.log().add(s2.log()).mul(-1 / 4)).exp());
	} else {
		throw new TypeError("Quat, quadrangle: Neplatný parametr: musí být Quat");
	}
};

/**
 * 
 * @param  {Quat} q1  kvaternion
 * @param  {Quat} q2  kvaternion
 * @param  {Quat} q3  kvaternion
 * @param  {number} t
 * @return {Quat} nová instance Quat
 * @throws {TypeError} If q1, q2, q3 nejsou Quat nebo t není číslo
 */
Quat.prototype.squad2 = function(q1, q2, q3, t) {
	if (arguments.length !== 4) {
		throw new TypeError("Quat, squad2: Neplatný počet parametrů: musí být 4");
	} else if (q1 instanceof Quat && q2 instanceof Quat && q3 instanceof Quat && typeof t === "number") {
		let s1 = this.quadrangle(q1, q2);
		let s2 = q2.quadrangle(this, q3);
		return new Quat(this.slerp(q2, t).slerp(s1.slerp(s2, t), (2 * t * (1 - t))));
	} else {
		throw new TypeError("Quat, squad2: Neplatný parametr: musí být Quat a číslo");
	}
};

/**
 * Výpis hodnot do konzole
 * @return {Quat} reference na volanou instanci pro možné další řetězení
 */
Quat.prototype.c = function() {
	console.log(this);
	return this;
};


/**
 * Objekt pro práci s maticemi 3×3
 * @param {Vec3D,Mat3,Mat4} v1
 * @param {Vec3D} v2
 * @param {Vec3D} v3
 * @constructor
 */
var Mat3 = function(v1, v2, v3) {
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
 * Vytváří jednotkovou matici 3×3
 * @augments {Mat3}
 */
var Mat3Identity = function() {
	this.mat = new ZeroArray(3);
	for (let i = 0; i < 3; i++) {
		this.mat[i][i] = 1;
	}
};
Mat3Identity.prototype = Object.create(Mat3.prototype);
Mat3Identity.prototype.constructor = Mat3Identity;
Mat3Identity.prototype.parent = Mat3;

/**
 * Vytváří transformační matici 3×3 pro rotaci kolem osy X ve 3D
 * @augments {Mat3}
 * @param {number} alpha úhel rotace v radiánech
 */
var Mat3RotX = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotX: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 3×3 pro rotaci kolem osy Y ve 3D
 * @augments {Mat3}
 * @param {number} alpha úhel rotace v radiánech
 */
var Mat3RotY = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotY: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 3×3 pro rotaci kolem osy Z ve 3D
 * @augments {Mat3}
 * @param {number} alpha úhel rotace v radiánech
 */
var Mat3RotZ = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat3RotZ: Neplatný parametr: musí být číslo");
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
 * Přičtení matice 3×3
 * @param {Mat3} m matice 3×3
 * @return {Mat3} nová instance Mat3
 * @throws {TypeError} If m není Mat3
 */
Mat3.prototype.add = function(m) {
	if (m instanceof Mat3) {
		let hlp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				hlp.mat[i][j] = this.mat[i][j] + m.mat[i][j];
			}
		}
		return hlp;
	} else {
		throw new TypeError("Mat3, add: Neplatný parametr: musí být Mat3");
	}
};

/**
 * Násobení matice skalárem, Přinásobení matice 3×3 zprava
 * @param  {number,Mat3} m skalár, matice 3×3
 * @return {Mat3} nová instance Mat3
 * @throws {TypeError} If m není číslo ani Mat3
 */
Mat3.prototype.mul = function(m) {
	if (typeof m === "number") {
		let hlp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				hlp.mat[i][j] = this.mat[i][j] * m;
			}
		}
		return hlp;
	} else if (m instanceof Mat3) {
		let hlp = new Mat3();
		for (let i = 0; i < 3; i++) {
			for (let j = 0; j < 3; j++) {
				let sum = 0;
				for (let k = 0; k < 3; k++) {
					sum += this.mat[i][k] * m.mat[k][j];
				}
				hlp.mat[i][j] = sum;
			}
		}
		return hlp;
	} else {
		throw new TypeError("Mat3, mul: Neplatný parametr: musí být číslo nebo Mat3");
	}
};

/**
 * Transponování matice 3×3
 * @return {Mat3} nová instance Mat3
 */
Mat3.prototype.transpose = function() {
	let hlp = new Mat3();
	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			hlp.mat[i][j] = this.mat[j][i];
		}
	}
	return hlp;
};

/**
 * Determinant matice 3×3
 * @return {number} determinant
 */
Mat3.prototype.det = function() {
	return this.mat[0][0]*(this.mat[1][1]*this.mat[2][2] - this.mat[2][1]*this.mat[1][2]) -
			this.mat[0][1]*(this.mat[1][0]*this.mat[2][2] - this.mat[2][0]*this.mat[1][2]) +
			this.mat[0][2]*(this.mat[1][0]*this.mat[2][1] - this.mat[2][0]*this.mat[1][1]);
};

/**
 * Inverzní matice 3×3
 * @return {Mat3} nová instance Mat3
 */
Mat3.prototype.inverse = function() {
	let hlp = new Mat3();
	let det = 1.0/this.det();
	hlp.mat[0][0] = det*(this.mat[1][1]*this.mat[2][2] - this.mat[1][2]*this.mat[2][1]);
	hlp.mat[0][1] = det*(this.mat[0][2]*this.mat[2][1] - this.mat[0][1]*this.mat[2][2]);
	hlp.mat[0][2] = det*(this.mat[0][1]*this.mat[1][2] - this.mat[0][2]*this.mat[1][1]);

	hlp.mat[1][0] = det*(this.mat[1][2]*this.mat[2][0] - this.mat[1][0]*this.mat[2][2]);
	hlp.mat[1][1] = det*(this.mat[0][0]*this.mat[2][2] - this.mat[0][2]*this.mat[2][0]);
	hlp.mat[1][2] = det*(this.mat[0][2]*this.mat[1][0] - this.mat[0][0]*this.mat[1][2]);

	hlp.mat[2][0] = det*(this.mat[1][0]*this.mat[2][1] - this.mat[1][1]*this.mat[2][0]);
	hlp.mat[2][1] = det*(this.mat[0][1]*this.mat[2][0] - this.mat[0][0]*this.mat[2][1]);
	hlp.mat[2][2] = det*(this.mat[0][0]*this.mat[1][1] - this.mat[0][1]*this.mat[1][0]);
	return hlp;
};

/**
 * Výpis matice do konzole
 * @return {Mat3} reference na volanou instanci pro možné další řetězení
 */
Mat3.prototype.c = function() {
	console.log("%cMat3", 'font-style: italic; font-weight: bold');
	for (let i = 0; i < 3; i++) {
		let x = "";
		for (let j = 0; j < 3; j++) {
			x += this.mat[i][j] + ", ";
		}
		console.log(x.substring(0, x.length-2));
	}
	return this;
};


/**
 * Objekt pro práci s maticemi 4×4
 * @param {Point3D,Mat4} p1
 * @param {Point3D} p2
 * @param {Point3D} p3
 * @param {Point3D} p4
 * @constructor
 */
var Mat4 = function(p1, p2, p3, p4) {
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
 * Vytváří jednotkou matici 4×4
 * @augments {Mat4}
 */
var Mat4Identity = function() {
	this.mat = new ZeroArray(4);
	for (let i = 0; i < 4; i++) {
		this.mat[i][i] = 1.0;
	}
};
Mat4Identity.prototype = Object.create(Mat4.prototype);
Mat4Identity.prototype.constructor = Mat4Identity;
Mat4Identity.prototype.parent = Mat4;

/**
 * Vytváří transformační matici 4×4 pro ortogonální deformaci zobrazovacího objemu
 * @augments {Mat4}
 * @param {number} w  šířka plátna
 * @param {number} h  výška plátna
 * @param {number} zn blízké z
 * @param {number} zf vzdálené z
 * @throws {TypeError} If některý z parametrů není zadán nebo není číslem
 */
var Mat4OrthoRH = function(w, h, zn, zf) {
	if (arguments.length !== 4) {
		throw new TypeError("Mat4OrthoRH: Neplatný počet parametrů: musí být 4");
	}
	if (typeof w !== "number" || typeof h !== "number" || typeof zn !== "number" || typeof zf !== "number") {
		throw new TypeError("Mat4OrthoRH: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro perspektivní deformaci zobrazovacího objemu
 * @augments {Mat4}
 * @param {number} alpha zorný úhel
 * @param {number} k     poměr šířka/výška plátna
 * @param {number} zn    blízké z
 * @param {number} zf    vzdálené z
 * @throws {TypeError} If některý z parametrů není zadán nebo není číslem
 */
var Mat4PerspRH = function(alpha, k, zn, zf) {
	if (arguments.length !== 4) {
		throw new TypeError("Mat4PerspRH: Neplatný počet parametrů: musí být 4");
	}
	if (typeof alpha !== "number" || typeof k !== "number" || typeof zn !== "number" || typeof zf !== "number") {
		throw new TypeError("Mat4PerspRH: Neplatný parametr: musí být číslo");
	}

	let h = (1.0 / Math.tan(alpha / 2.0));
	let w = k * h;
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
 * Vytváří transformační matici 4×4 pro rotaci kolem osy X ve 3D
 * @augments {Mat4}
 * @param {number} alpha úhel rotace v radiánech
 * @throws {TypeError} If alpha není číslo
 */
var Mat4RotX = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotX: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro rotaci kolem osy Y ve 3D
 * @augments {Mat4}
 * @param {number} alpha úhel rotace v radiánech
 * @throws {TypeError} If alpha není číslo
 */
var Mat4RotY = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotY: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro rotaci kolem osy Z ve 3D
 * @augments {Mat4}
 * @param {number} alpha úhel rotace v radiánech
 * @throws {TypeError} If alpha není číslo
 */
var Mat4RotZ = function(alpha) {
	if (typeof alpha !== "number") {
		throw new TypeError("Mat4RotZ: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro rotaci kolem os X, Y, Z ve 3D
 * @augments {Mat4}
 * @param {number} alpha úhel rotace v radiánech
 * @param {number} beta  úhel rotace v radiánech
 * @param {number} gama  úhel rotace v radiánech
 * @throws {TypeError} If některý z argumentů není zadán nebo není číslem
 */
var Mat4RotXYZ = function(alpha, beta, gama) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4RotXYZ: Neplatný počet parametrů: musí být 3");
	}
	if (typeof alpha !== "number" || typeof beta !== "number" || typeof gama !== "number") {
		throw new TypeError("Mat4RotXYZ: Neplatný parametr: musí být číslo");
	}

	this.mat = new Mat4RotX(alpha).mul(new Mat4RotY(beta)).mul(new Mat4RotZ(gama)).mat;
};
Mat4RotXYZ.prototype = Object.create(Mat4.prototype);
Mat4RotXYZ.prototype.constructor = Mat4RotXYZ;
Mat4RotXYZ.prototype.parent = Mat4;

/**
 * Vytváří transformační matici 4×4 pro změnu měřítka ve 3D
 * @augments {Mat4}
 * @param {number} x zvětšení/zmenšení na ose x
 * @param {number} y zvětšení/zmenšení na ose y
 * @param {number} z zvětšení/zmenšení na ose z
 * @throws {TypeError} If některý z argumentů není zadán nebo není číslem
 */
var Mat4Scale = function(x, y, z) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4Scale: Neplatný počet parametrů: musí být 3");
	}
	if (typeof x !== "number" || typeof y !== "number" || typeof z !== "number") {
		throw new TypeError("Mat4Scale: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro translaci ve 3D
 * @augments {Mat4}
 * @param {number} x posunutí na ose x
 * @param {number} y posunutí na ose y
 * @param {number} z posunutí na ose z
 * @throws {TypeError} If některý z argumentů není zadán nebo není číslem
 */
var Mat4Transl = function(x, y, z) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4Transl: Neplatný počet parametrů: musí být 3");
	}
	if (typeof x !== "number" || typeof y !== "number" || typeof z !== "number") {
		throw new TypeError("Mat4Transl: Neplatný parametr: musí být číslo");
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
 * Vytváří transformační matici 4×4 pro pohledovou transformaci ve 3D
 * @augments {Mat4}
 * @param {number} e vektor pozice pozorovatele
 * @param {number} v vektor pohledu
 * @param {number} u up vektor
 * @throws {TypeError} If některý z argumentů není zadán nebo není číslem
 */
var Mat4ViewRH = function(e, v, u) {
	if (arguments.length !== 3) {
		throw new TypeError("Mat4ViewRH: Neplatný počet parametrů: musí být 3");
	}
	if (!(e instanceof Vec3D) || !(v instanceof Vec3D) || !(u instanceof Vec3D)) {
		throw new TypeError("Mat4ViewRH: Neplatný parametr: musí být Vec3D");
	}

	this.mat = new Mat4Identity().mat;
	let z = v.mul(-1.0).normalized();
	let x = u.cross(z).normalized();
	let y = z.cross(x);

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
 * Sčítání matic 4×4
 * @param {Mat4} m matice 4×4
 * @return {Mat4} nová instance Mat4
 * @throws {TypeError} If m není Mat4
 */
Mat4.prototype.add = function(m) {
	if (m instanceof Mat4) {
		let hlp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				hlp.mat[i][j] = this.mat[i][j] + m.mat[i][j];
			}
		}
		return hlp;
	} else {
		throw new TypeError("Mat4, add: Neplatný parametr: musí být Mat4");
	}
};

/**
 * Násobení matice 4×4 skalárem, Násobení maticí 4×4 zprava
 * @param  {number,Mat4} m skalár, matice 4×4
 * @return {Mat4} nová instance Mat4
 * @throws {TypeError} If m není číslo ani Mat4
 */
Mat4.prototype.mul = function(m) {
	if (typeof m === "number") {
		let tmp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				tmp.mat[i][j] = this.mat[i][j] * m;
			}
		}
		return tmp;
	} else if (m instanceof Mat4) {
		let tmp = new Mat4();
		for (let i = 0; i < 4; i++) {
			for (let j = 0; j < 4; j++) {
				let sum = 0;
				for (let k = 0; k < 4; k++) {
					sum += this.mat[i][k] * m.mat[k][j];
				}
				tmp.mat[i][j] = sum;
			}
		}
		return tmp;
	} else {
		throw new TypeError("Mat4, mul: Neplatný parametr: musí být číslo nebo Mat4");
	}
};

/**
 * Transponování matice 4×4
 * @return {Mat4} nová instance Mat4
 */
Mat4.prototype.transpose = function() {
	let tmp = new Mat4();
	for (let i = 0; i < 4; i++) {
		for (let j = 0; j < 4; j++) {
			tmp.mat[i][j] = this.mat[j][i];
		}
	}
	return tmp;
};

/**
 * Převod Mat4 na Mat3, odstranění posledního sloupce a řádku
 * @return {Mat3} nová instance Mat3
 */
Mat4.prototype.toMat3 = function() {
	let tmp = new Mat3();
	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			tmp.mat[i][j] = this.mat[i][j];
		}
	}
	return tmp;
};

/**
 * Výpis matice do konzole
 * @return {Mat4} reference na volanou instanci pro možné další řetězení
 */
Mat4.prototype.c = function() {
	console.log("%cMat4", 'font-style: italic; font-weight: bold');
	for (let i = 0; i < 4; i++) {
		let x = "";
		for (let j = 0; j < 4; j++) {
			x += this.mat[i][j] + ", ";
		}
		console.log(x.substring(0, x.length-2));
	}
	return this;
};


/**
 * Objekt pro nastavení pohledové transformace
 * @constructor
 */
var Camera = function() {
	this.azimuth = 0.0;
	this.zenith = 0.0;
	this.radius = 1.0;
	this.xy = true;
	this.pos = new Vec3D(0.0, 0.0, 0.0);
	this.firstPerson = true; // true -> 1. osoba, false -> 3. osoba

	this.eye = new Vec3D();
	this.eyeVector = new Vec3D();
	this.view = new Mat4();

	this.computeMatrix();
};

/**
 * Přepočítání eye, eyeVector a view
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
 * Přičtení azimutu
 * @param {number} ang
 */
Camera.prototype.addAzimuth = function(ang) {
	if (typeof ang === "number") {
		this.azimuth += ang;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, addAzimuth: Neplatný parametr: musí být číslo");
	}
};

/**
 * Nastavení azimutu
 * @param {number} ang
 * @throws {TypeError} If ang není číslo
 */
Camera.prototype.setAzimuth = function(ang) {
	if (typeof ang === "number") {
		this.azimuth = ang;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setAzimuth: Neplatný parametr: musí být číslo");
	}
};

/**
 * Přičtení zenitu
 * @param {number} ang
 * @throws {TypeError} If ang není číslo
 */
Camera.prototype.addZenith = function(ang) {
	if (typeof ang === "number") {
		if (Math.abs(this.zenith + ang) <= Math.PI / 2) {
			this.zenith += ang;
			this.computeMatrix();
		}
	} else {
		throw new TypeError("Camera, addZenith: Neplatný parametr: musí být číslo");
	}
};

/**
 * Nastavení zenitu
 * @param {number} ang
 * @throws {TypeError} If ang není číslo
 */
Camera.prototype.setZenith = function(ang) {
	if (typeof ang === "number") {
		this.zenith = ang;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setZenith: Neplatný parametr: musí být číslo");
	}
};

/**
 * Přičení radiusu
 * @param {number} dist
 * @throws {TypeError} If dist není číslo
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
		throw new TypeError("Camera, addRadius: Neplatný parametr: musí být číslo");
	}
};

/**
 * Nastavení radiusu
 * @param {number} dist
 * @throws {TypeError} If dist není číslo
 */
Camera.prototype.setRadius = function(dist) {
	if (typeof dist === "number") {
		this.radius = dist;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setRadius: Neplatný parametr: musí být číslo");
	}
};

/**
 * Vynásobit radius
 * @param  {number} scale
 * @throws {TypeError} If scale není číslo
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
		throw new TypeError("Camera, mulRadius: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb dopředu
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.forward = function(speed) {
	if (typeof speed === "number") {
		if (!this.xy) {
			this.pos = this.pos.add(new Vec3D(
				(Math.cos(this.azimuth - Math.PI / 2) * Math.cos(this.zenith + Math.PI)),
				(Math.sin(this.azimuth - Math.PI / 2) * Math.cos(this.zenith + Math.PI)),
				Math.sin(this.zenith)).mul(speed)
			);
		} else {
			this.pos = this.pos.add(
				new Vec3D(
					Math.cos(this.azimuth - Math.PI / 2),
					Math.sin(this.azimuth - Math.PI / 2),
					0.0
				).mul(-speed)
			);
		}
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, forward: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb dozadu
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.backward = function(speed) {
	if (typeof speed === "number") {
		this.forward(-speed);
	} else {
		throw new TypeError("Camera, backward: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb doprava
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.right = function(speed) {
	if (typeof speed === "number") {
		this.pos = this.pos.add(
			new Vec3D(Math.cos(this.azimuth), Math.sin(this.azimuth), 0).mul(speed)
		);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, right: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb doleva
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.left = function(speed) {
	if (typeof speed === "number") {
		this.right(-speed);
	} else {
		throw new TypeError("Camera, left: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb dolů
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.down = function(speed) {
	if (typeof speed === "number") {
		this.pos.z -= speed;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, down: Neplatný parametr: musí být číslo");
	}
};

/**
 * Pohyb nahoru
 * @param  {number} speed
 * @throws {TypeError} If speed není číslo
 */
Camera.prototype.up = function(speed) {
	if (typeof speed === "number") {
		this.pos.z += speed;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, up: Neplatný parametr: musí být číslo");
	}
};

/**
 * Přičtení pozice
 * @param  {Vec3D} dir
 * @throws {TypeError} If dir není Vec3D
 */
Camera.prototype.move = function(dir) {
	if (dir instanceof Vec3D) {
		this.pos = this.pos.add(dir);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, move: Neplatný parametr: musí být Vec3D");
	}
};

/**
 * Nastavení pozice
 * @param  {Vec3D} apos
 * @throws {TypeError} If apos není Vec3D
 */
Camera.prototype.setPosition = function(apos) {
	if (apos instanceof Vec3D) {
		this.pos = new Vec3D(apos);
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setPosition: Neplatný parametr: musí být Vec3D");
	}
};

/**
 * Nastavení první osoby
 * @param  {boolean} fp true-> 1. osoba; false-> 3. osoba
 * @throws {TypeError} If fp není boolean
 */
Camera.prototype.setFirstPerson = function(fp) {
	if (typeof fp === "boolean") {
		this.firstPerson = fp;
		this.computeMatrix();
	} else {
		throw new TypeError("Camera, setFirstPerson: Neplatný parametr: musí být boolean");
	}
};

/**
 * Objekt pro práci s barvami
 * @param {number,Col,Point3D} ar červená složka, barva, Point3D
 * @param {number} ag zelená složka
 * @param {number} ab modrá složka
 * @param {number} aa alfa složka
 * @constructor
 * @throws {TypeError} If parametry nesplňují ani jednu podmínku
 */
var Col = function (ar, ag, ab, aa) {
	// založit barvu na jiné barvě
	if (ar instanceof Col) {
		this.r = ar.r;
		this.g = ar.g;
		this.b = ar.b;
		this.a = ar.a;
	//založit barvu na bodu
	} else if (ar instanceof Point3D) {
		this.r = ar.x;
		this.g = ar.y;
		this.b = ar.z;
		this.a = ar.w;
	// jsou 3 nebo 4 parametry a jsou to celá čísla, pak se předpokládá interval <0; 255>
	// !!! pozor např. r=1, g=0, b=1 projde podmínkami
	} else if ((arguments.length === 3 || arguments.length === 4) &&
			this.isInt(ar) && this.isInt(ag) && this.isInt(ab)) {
		this.r = ar / 255.0;
		this.g = ag / 255.0;
		this.b = ab / 255.0;
		this.a = (this.isInt(aa)) ? aa / 255.0 : 1.0;
	// jsou 3 nebo 4 parametry a nejsou to celá čísla, pak se předpokládá interval <0; 1>
	} else if (arguments.length === 3 || arguments.length === 4) {
		this.r = (typeof ar === "number") ? ar : 1.0;
		this.g = (typeof ag === "number") ? ag : 1.0;
		this.b = (typeof ab === "number") ? ab : 1.0;
		this.a = (typeof aa !== "undefined" && typeof aa === "number") ? aa : 1.0;
	} else {
		throw new TypeError("Col: Neplatný parametr.");
	}
};

/**
 * Je číslo celé (integer)
 * @param  {number}  n číslo
 * @return {Boolean}   true, pokud celé číslo, jinak false
 */
Col.prototype.isInt = function(n) {
	return (n % 1 === 0);
};

/**
 * Přičíst barvu
 * @param  {Col} c barva pro přičtení
 * @return {Col}   nová instance Col
 * @throws {TypeError} If c není Col
 */
Col.prototype.addna = function(c) {
	if (c instanceof Col) {
		return new Col(this.r + c.r, this.g + c.g, this.b + c.b);
	} else {
		throw new TypeError("Col, addna: Neplatný parametr: musí být Col");
	}
};

/**
 * Vynásobit číslem
 * @param  {number} x číslo
 * @return {Col}      nová instance Col
 * @throws {TypeError} If x není číslo
 */
Col.prototype.mulna = function(x) {
	if (typeof x === "number") {
		return new Col(this.r * x, this.g * x, this.b * x);
	} else {
		throw new TypeError("Col, mulna: Neplatný parametr: musí být číslo");
	}
};

/**
 * Přičíst barvu včetně alfa složky
 * @param {Col} c barva pro přičtení
 * @return {Col}  nová instance Col
 * @throws {TypeError} If c není Col
 */
Col.prototype.add = function(c) {
	if (c instanceof Col) {
		return new Col(this.r + c.r, this.g + c.g, this.b + c.b, this.a + c.a);
	} else {
		throw new TypeError("Col, add: Neplatný parametr: musí být Col");
	}
};

/**
 * Vynásobit barvu včetně alfa složky
 * @param  {number,Col} c číslo nebo Col
 * @return {Col}          nová instance Col
 * @throws {TypeError} If c není Col ani číslo
 */
Col.prototype.mul = function(c) {
	if (c instanceof Col) {
		return new Col(this.r * c.r, this.g * c.g, this.b * c.b, this.a * c.a);
	} else if (typeof c === "number") {
		return new Col(this.r * c, this.g * c, this.b * c, this.a * c);
	} else {
		throw new TypeError("Col, mul: Neplatný parametr: musí být Col nebo číslo");
	}
};

/**
 * Gamma korekce
 * @param  {number} gamma
 * @return {Col}          nová instance Col
 * @throws {TypeError} If gamma není číslo
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
		throw new TypeError("Col, gamma: Neplatný parametr: musí být číslo");
	}
};

/**
 * Saturace barvy
 * @return {Col} nová instance Col
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
 * Získání barvy jako 1 číslo, bez afla složky
 * @return {number} barva
 */
Col.prototype.getRGB = function() {
	return ((this.r * 255.0) << 16) | ((this.g * 255.0) << 8) | (this.b * 255.0);
};

/**
 * Získání barvy jako 1 číslo, s alfa složkou
 * @return {number} barva
 */
Col.prototype.getARGB = function() {
	return ((this.a * 255.0) << 24) | ((this.r * 255.0) << 16) | ((this.g * 255.0) << 8) | (this.b * 255.0);
};

/**
 * Výpis barvy do konzole
 * @return {Col} reference na volanou instanci pro možné další řetězení
 */
Col.prototype.c = function() {
	console.log(this);
	return this;
};

/**
 * Objekt pro práci s kubikami - fergusonovy, beziérovy a coonsovy
 * @param {number} type 1-> ferguson, 2-> coons, ostatní->beziér
 * @constructor
 */
var Kubika = function(type) {
	//bazova matice
	this.bm = new Mat4();
	//matice ridicich bodu
	this.rb;

	switch (typeof type === "number" ? type : 0) {
		case 1:// ferguson
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

		case 2:// coons
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

		case 0:// bezier
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
 * Inicializace pomocí zadané čtveřice řídících bodů
 * @param  {Point3D} b1 řídící bod
 * @param  {Point3D} b2 řídící bod
 * @param  {Point3D} b3 řídící bod
 * @param  {Point3D} b4 řídící bod
 * @throws {TypeError} If některý z paramterů není zadán nebo není Point3D
 */
Kubika.prototype.init = function(b1, b2, b3, b4) {
	if (arguments.length !== 4) {
		throw new TypeError("Kubika, init: Neplatný počet parametrů: musí být 4");
	} else if (b1 instanceof Point3D && b2 instanceof Point3D && b3 instanceof Point3D && b4 instanceof Point3D) {
		if (this.type === 1) {
			this.rb = new Mat4(b1, b4, b2.sub(b1), b4.sub(b3));
		} else {
			this.rb = new Mat4(b1, b2, b3, b4);
		}
		this.rb = this.bm.mul(this.rb);
	} else {
		throw new TypeError("Kubika, init: Neplatný parametr: musí být Point3D");
	}
};

/**
 * Výpočet souřadnice bodu kubiky podle parametru t
 * @param  {number} t parametr z intervalu <0; 1>; pokud mimo, tak oříznuto
 * @return {Point3D}  výsledný bod
 * @throws {TypeError} If t není číslo
 */
Kubika.prototype.compute = function(t) {
	if (typeof t === "number") {
		if (t > 1) t = 1;
		if (t < 0) t = 0;

		let res = new Point3D(t * t * t, t * t, t, 1);

		res = res.mul(this.rb);
		res.w = 1;
		return res;
	} else {
		throw new TypeError("Kubika, compute: Neplatný parametr: musí být číslo");
	}
};

/**
 * Objekt pro práci s bikubikami
 * @param {number} type 1-> ferguson, 2-> coons, ostatní->beziér
 * @constructor
 */
var Bikubika = function(type) {
	//Point3D
	this.u1; this.u2; this.u3; this.u4;

	this.k1 = new Kubika(type);
	this.k2 = new Kubika(type);
	this.k3 = new Kubika(type);
	this.k4 = new Kubika(type);
	this.k5 = new Kubika(type);
};

/**
 * Inicializace pomocí 4×4 řídících bodů
 * @param  {Point3D} b11 řídící bod
 * @param  {Point3D} b12 řídící bod
 * @param  {Point3D} b13 řídící bod
 * @param  {Point3D} b14 řídící bod
 * @param  {Point3D} b21 řídící bod
 * @param  {Point3D} b22 řídící bod
 * @param  {Point3D} b23 řídící bod
 * @param  {Point3D} b24 řídící bod
 * @param  {Point3D} b31 řídící bod
 * @param  {Point3D} b32 řídící bod
 * @param  {Point3D} b33 řídící bod
 * @param  {Point3D} b34 řídící bod
 * @param  {Point3D} b41 řídící bod
 * @param  {Point3D} b42 řídící bod
 * @param  {Point3D} b43 řídící bod
 * @param  {Point3D} b44 řídící bod
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
 * Výpočet souřadnice bodu bikubiky podle paramterů u, v
 * @param  {number} u parametr z intervalu <0; 1>; pokud mimo, tak oříznuto
 * @param  {number} v parametr z intervalu <0; 1>; pokud mimo, tak oříznuto
 * @return {Point3D}  souřadnice bodu
 * @throws {TypeError} If u nebo v není číslo
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
		throw new TypeError("Bikubika, compute: Neplatný parametr: musí být číslo");
	}
};
