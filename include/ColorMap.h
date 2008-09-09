#ifndef COLORMAP_H_
#define COLORMAP_H_

/**
 * This class has been obtained from 
 * http://www.koders.com/cpp/fidCFEC52F6E8D5CAF77FBCDB42FDB45A69EC48677B.aspx?s=colorMapFunc#L153
 * under the GNU licence.
 * Minor modifications by:
 * @author: S.Daraszewicz
 * @date: 28.08.08
 */

typedef void (*colorMapFunc)(unsigned int*,float,float,float);

class ColorMap{
 public:
  static void colorMap(unsigned int *rgb,float value,float min,float max);
  static void hotColorMap(unsigned int *rgb,float value,float min,float max);
  static void coldColorMap(unsigned int *rgb,float value,float min,float max);
  static void jetColorMap(unsigned int *rgb,float value,float min,float max);
  static void cyclicColorMap(unsigned int *rgb,float value,float min,float max);
  static void randColorMap(unsigned int *rgb,float value,float min,float max);
  static void grayColorMap(unsigned int *rgb,float value,float min,float max);
  static void blueColorMap(unsigned int *rgb,float value,float min,float max);
  static colorMapFunc selectColorMap(int cmp);
  static int RGB2HEX(int red, int green, int blue);
};

#endif
