#ifndef MINMAX_H
#define MINMAX_H 1

class dMinMax {
public:
    double min;
    double max;
    void set(double v){
	if(v<min)
	    min=v;
	if(v>max)
	    max=v;
    }
};

class iMinMax {
public:
    int min;
    int max;
    void set(int v){
	if(v<min)
	    min=v;
	if(v>max)
	    max=v;
    }
};

#endif
