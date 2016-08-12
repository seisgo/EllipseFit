#include "DirectEllipseFit.h"
#include <math.h>

Ellipse::Ellipse()
{
    algeFlag = false;
    a = b = c = d = e = f = 0;

    geomFlag = false;
    cx = cy = 0;
    rl = rs = 0;
    phi = 0;
}

void Ellipse::alge2geom()
{
    if(!algeFlag)
        return;

    double tmp1 = b*b - 4*a*c;
    double tmp2 = sqrt((a-c)*(a-c)+b*b);
    double tmp3 = a*e*e + c*d*d - b*d*e + tmp1*f;

    double r1 = -sqrt(2*tmp3*(a+c+tmp2)) / tmp1;
    double r2 = -sqrt(2*tmp3*(a+c-tmp2)) / tmp1;
    rl = r1>=r2 ? r1 : r2;
    rs = r1<=r2 ? r1 : r2;

    cx = (2*c*d - b*e) / tmp1;
    cy = (2*a*e - b*d) / tmp1;

    phi = 0.5 * atan2(b, a-c);
    if(r1>r2)
        phi += M_PI_2;

    geomFlag = true;
}

void Ellipse::geom2alge()
{
    if(!geomFlag)
        return;

    a = rl*rl*sin(phi)*sin(phi) + rs*rs*cos(phi)*cos(phi);
    b = 2 * (rs*rs - rl*rl) * sin(phi) * cos(phi);
    c = rl*rl*cos(phi)*cos(phi) + rs*rs*sin(phi)*sin(phi);
    d = -2*a*cx - b*cy;
    e = -b*cx - 2*c*cy;
    f = a*cx*cx + b*cx*cy + c*cy*cy - rl*rl*rs*rs;

    algeFlag = true;
}
