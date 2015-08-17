

typedef struct complex{
	double real;
	double imag;
}complex;

struct complex add(struct complex c1, struct complex c2){
	double realSum = c1.real + c2.real;
	double imagSum = c1.imag + c2.imag;

	complex result;
	result.real = realSum;
	result.imag = imagSum;
	return result;
}

struct complex subtract(struct complex c1, struct complex c2){
	double realDiff = c1.real - c2.real;
	double imagDiff = c1.imag - c2.imag;

	complex result;
	result.real = realDiff;
	result.imag = imagDiff;
	return result;
}

struct complex mulitply(struct complex c1, struct complex c2){
	double aSquared = c1.real * c2.real;
	double twoAB = (c1.real * c2.imag) + (c1.imag * c2.real);
	double bSquared = -(c1.imag * c2.imag);

	complex result;
	result.real = aSquared + bSquared;
	result.imag = twoAB;
	return result;
}

double squaredNorm(struct complex c){
	double result = (c.real * c.real) + (c.imag * c.imag);
	return result;
}

