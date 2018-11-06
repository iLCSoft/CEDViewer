#include <math.h>
#include <stdlib.h>
#include "ColorMap.h"

/**
 * This class has been obtained from 
 * http://www.koders.com/cpp/fidCFEC52F6E8D5CAF77FBCDB42FDB45A69EC48677B.aspx?s=colorMapFunc#L153
 * under the GNU licence.
 * Minor modifications by:
 * @author: S.Daraszewicz
 * @date: 28.08.08
 */

void ColorMap::jetColorMap(unsigned int *rgb,float value,float min,float max)
{
  unsigned int c1=144;
  float max4=(max-min)/4;
  value-=min;
  if(value==HUGE_VAL)
    {rgb[0]=rgb[1]=rgb[2]=255;}
  else if(value<0)
    {rgb[0]=rgb[1]=rgb[2]=0;}
  else if(value<max4)
    {rgb[0]=0;rgb[1]=0;rgb[2]=c1+(unsigned int)((255-c1)*value/max4);}
  else if(value<2*max4)
    {rgb[0]=0;rgb[1]=(unsigned int)(255*(value-max4)/max4);rgb[2]=255;}
  else if(value<3*max4)
    {rgb[0]=(unsigned int)(255*(value-2*max4)/max4);rgb[1]=255;rgb[2]=255-rgb[0];}
  else if(value<max)
    {rgb[0]=255;rgb[1]=(unsigned int)(255-255*(value-3*max4)/max4);rgb[2]=0;}
  else {rgb[0]=255;rgb[1]=rgb[2]=0;}
}

void ColorMap::hotColorMap(unsigned int *rgb,float value,float min,float max)
{
  float max3=(max-min)/3;
  value-=min;
  if(value==HUGE_VAL)
    {rgb[0]=rgb[1]=rgb[2]=255;}
  else if(value<0)
    {rgb[0]=rgb[1]=rgb[2]=0;}
  else if(value<max3)
    {rgb[0]=(unsigned int)(255*value/max3);rgb[1]=0;rgb[2]=0;}
  else if(value<2*max3)
    {rgb[0]=255;rgb[1]=(unsigned int)(255*(value-max3)/max3);rgb[2]=0;}
  else if(value<max)
    {rgb[0]=255;rgb[1]=255;rgb[2]=(unsigned int)(255*(value-2*max3)/max3);}
  else {rgb[0]=rgb[1]=rgb[2]=255;}
}

void ColorMap::coldColorMap(unsigned int *rgb,float value,float min,float max)
{
  float max3=(max-min)/3;
  value-=min;
  if(value==HUGE_VAL)
    {rgb[0]=rgb[1]=rgb[2]=255;}
  else if(value<0)
    {rgb[0]=rgb[1]=rgb[2]=0;}
  else if(value<max3)
    {rgb[0]=0;rgb[1]=0;rgb[2]=(unsigned int)(255*value/max3);}
  else if(value<2*max3)
    {rgb[0]=0;rgb[1]=(unsigned int)(255*(value-max3)/max3);rgb[2]=255;}
  else if(value<max)
    {rgb[0]=(unsigned int)(255*(value-2*max3)/max3);rgb[1]=255;rgb[2]=255;}
  else {rgb[0]=rgb[1]=rgb[2]=255;}
}

void ColorMap::blueColorMap(unsigned int *rgb,float value,float min,float max)
{
  value-=min;
  if(value==HUGE_VAL)
    {rgb[0]=rgb[1]=rgb[2]=255;}
  else if(value<0)
    {rgb[0]=rgb[1]=rgb[2]=0;}
  else if(value<max)
    {rgb[0]=0;rgb[1]=0;rgb[2]=(unsigned int)(255*value/max);}
  else {rgb[0]=rgb[1]=0;rgb[2]=255;}
}

void positiveColorMap(unsigned int *rgb,float value,float min,float max)
{
  value-=min;
  max-=min;
  value/=max;

  if(value<0){
  rgb[0]=rgb[1]=rgb[2]=0;
    return;
  }
  if(value>1){
  rgb[0]=rgb[1]=rgb[2]=255;
  return;
  }

  rgb[0]=192;rgb[1]=0;rgb[2]=0;
  rgb[0]+=(unsigned int)(63*value);
  rgb[1]+=(unsigned int)(255*value);
  if(value>0.5)
  rgb[2]+=(unsigned int)(255*2*(value-0.5));
}

void negativeColorMap(unsigned int *rgb,float value,float min,float max)
{
  value-=min;
  max-=min;
  value/=max;

  rgb[0]=0;rgb[1]=0;rgb[2]=0;
  if(value<0) return;
  if(value>1){
  rgb[1]=rgb[2]=255;
  return;
  }

  rgb[1]+=(unsigned int)(255*value);
  if(value>0.5)
  rgb[2]+=(unsigned int)(255*2*(value-0.5));

}

void ColorMap::colorMap(unsigned int *rgb,float value,float min,float max)
{
  if(value>0) 
    positiveColorMap(rgb,value,0,max);
  else 
    negativeColorMap(rgb,value,min,0);
/*
  if(value>0) 
    hotColorMap(rgb,value,min,max);
  else 
    coldColorMap(rgb,value,min,max);
	*/
}

void ColorMap::cyclicColorMap(unsigned int *rgb,float value,float min,float max)
{
  float max3=(max-min)/3;
  value-=(max-min)*(float)floor((value-min)/(max-min));
  if(value<max3)
    {rgb[0]=(unsigned int)(255-255*value/max3);rgb[1]=0;rgb[2]=255-rgb[0];}
  else if(value<2*max3)
    {rgb[0]=0;rgb[1]=(unsigned int)(255*(value-max3)/max3);rgb[2]=255-rgb[1];}
  else if(value<max)
    {rgb[0]=(unsigned int)(255*(value-2*max3)/max3);rgb[1]=255-rgb[0];rgb[2]=0;}

}
void ColorMap::randColorMap(unsigned int *rgb,float value,float min,float max)
{
  srand((int)(65000*(value-min)/(max-min)));
  rgb[0]=(unsigned int)(255*rand());
  rgb[1]=(unsigned int)(255*rand());
  rgb[2]=(unsigned int)(255*rand());
}

void ColorMap::grayColorMap(unsigned int *rgb,float value,float min,float max)
{
  max-=min;
  value-=min;
  rgb[0]=rgb[1]=rgb[2]=(unsigned int)(255*value/max);
}

	
int ColorMap::RGB2HEX(int red, int green, int blue){
		
	int rgb = 0x000000; //default
    rgb = (red<<16) + (green<<8) + (blue); //bit shifting
	return rgb;
}

// NOT TESTED! ALPHA CHANNEL
//int ColorMap::RGBA2HEX(int red, int green, int blue, int alpha){
//		
//	int rgb = 0x000000; //default
//    rgb = (alpha<<24) + (red<<16) + (green<<8) + (blue); //bit shifting
//	return rgb;
//}

colorMapFunc ColorMap::selectColorMap(int cmp)
{
  int max=7;
  cmp=abs(cmp)%max;
  switch(cmp){
  case 0:
    return colorMap;
    break;
  case 1:
    return hotColorMap;
    break;
  case 2:
    return coldColorMap;
    break;
  case 3:
    return jetColorMap;
    break;
  case 4:
    return cyclicColorMap;
    break;
  case 5:
    return grayColorMap;
    break;
  case 6:
    return blueColorMap;
    break;
  default:
    return colorMap;
  };
}
// taken from https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both#comment66194776_14733008
RgbColor ColorMap::HsvToRgb(HsvColor in){
    double      hh, p, q, t, ff;
    long        i;
    RgbColor    out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = in.v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.v;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = in.v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.v;
        break;
    case 5:
    default:
        out.r = in.v;
        out.g = p;
        out.b = q;
        break;
    }
    return out;
}

// Convert a number between min in max to a color between blue(min) and red(max)
unsigned long ColorMap::NumberToTemperature(double value, double min, double max, double s, double v){
    if( value > max ){
        value = max;
    }else if ( value < min ){
        value = min;
    }

    HsvColor hsv;
    hsv.h = 240 - ( value - min ) / ( max - min ) * 240;
    hsv.s=s;
    hsv.v=v;

    RgbColor rgb = ColorMap::HsvToRgb(hsv);

    return ColorMap::RGB2HEX((int)(rgb.r*255),(int)(rgb.g*255),(int)(rgb.b*255));
}
