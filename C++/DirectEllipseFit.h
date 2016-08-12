/*******************************************************************************
 * ProName: DirectEllipseFit.h
 * Author:  Zhenyu Yuan
 * Data:    2016/7/27
 * -----------------------------------------------------------------------------
 * INSTRUCTION: Perform ellipse fitting by direct method
 * DEPENDANCE:  clapack and Qt
 * REFERENCE:
 *      (1) Fitzgibbon, A., et al. (1999). "Direct least square fitting of ellipses."
 *          IEEE Transactions on pattern analysis and machine intelligence 21(5):
 *          476-480.
 *      (2) http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/
 * CONVENTION:
 *      (1) Matrix expressed as QVector<QVector<T> >, fast direction is rows,
 *          slow direction is columns, that is, internal QVector indicates column
 *          vector of T, external QVector indicates row vector of column vectors
 *      (2) Matrix expressed as 1-order array, arranged in rows, that is, row by
 *          row.
 ******************************************************************************/
#ifndef DIRECTELLIPSEFIT_H
#define DIRECTELLIPSEFIT_H

#include <algo_global.h>
#include <QVector>
#include <clapack.h>
#include <QDebug>

class EP_ALGO_EXPORT Ellipse
{
public:
    Ellipse();
    /**
     * @brief alge2geom:    algebraic parameters to geometric parameters
     * @ref:    https://en.wikipedia.org/wiki/Ellipse#In_analytic_geometry
     *          http:homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/
     * @note:   The calculation of phi refer to wikipedia is not correct,
     *          refer to Bob Fisher's matlab program.
     *          What's more, the calculated geometric parameters can't back to
     *          initial algebraic parameters from geom2alge();
     */
    void alge2geom();
    /**
     * @brief geom2alge:    geometric parameters to algebraic parameters
     * @ref:    https://en.wikipedia.org/wiki/Ellipse#In_analytic_geometry
     */
    void geom2alge();

public:
    //algebraic parameters as coefficients of conic section
    float a, b, c, d, e, f;
    bool algeFlag;

    //geometric parameters
    float cx;   //centor in x coordinate
    float cy;   //centor in y coordinate
    float rl;   //semimajor: large radius
    float rs;   //semiminor: small radius
    float phi;  //azimuth angel in radian unit
    bool geomFlag;
};

template <typename T>
class DirectEllipseFit
{
public:
    DirectEllipseFit(const QVector<T> &xData, const QVector<T> &yData);
    Ellipse doEllipseFit();

private:
    T getMeanValue(const QVector<T> &data);
    T getMaxValue(const QVector<T> &data);
    T getMinValue(const QVector<T> &data);
    T getScaleValue(const QVector<T> &data);
    QVector<T> symmetricNormalize(const QVector<T> &data);
    //Make sure xData and yData are of same size
    QVector<T> dotMultiply(const QVector<T> &xData, const QVector<T> &yData);
    //Get n*6 design matrix D, make sure xData and yData are of same size
    QVector<QVector<T> > getDesignMatrix(const QVector<T> &xData,
                                         const QVector<T> &yData);
    //Get 6*6 constraint matrix C
    QVector<QVector<T> > getConstraintMatrix();
    //Get 6*6 scatter matrix S from design matrix
    QVector<QVector<T> > getScatterMatrix(const QVector<QVector<T> > &dMtrx);
    //Transpose matrix
    QVector<QVector<T> > transposeMatrix(const QVector<QVector<T> > &mtrx);
    //Do matrix multiplication, mtrx1: j*l; mtrx2: l*i; return: j*i
    QVector<QVector<T> > doMtrxMul(const QVector<QVector<T> > &mtrx1,
                                   const QVector<QVector<T> > &mtrx2);

    /**
     * @brief solveGeneralEigens:   Solve generalized eigensystem
     * @note        For real eiginsystem solving.
     * @param sMtrx:    6*6 square matrix in this application
     * @param cMtrx:    6*6 square matrix in this application
     * @param eigVV:    eigenvalues and eigenvectors, 6*7 matrix
     * @return  success or failure status
     */
    bool solveGeneralEigens(const QVector<QVector<T> > &sMtrx,
                            const QVector<QVector<T> > &cMtrx,
                            QVector<QVector<T> > &eigVV);
    //Convert matrix expression from nested QVector to 1-order array
    double *mtrx2array(const QVector<QVector<T> > &mtrx);

    /**
     * @brief calcEllipsePara:  calculate ellipse parameter form eigen information
     * @param eigVV:    eigenvalues and eigenvectors
     * @return ellipse parameter
     */
    Ellipse calcEllipsePara(const QVector<QVector<T> > &eigVV);

private:
    QVector<T> m_xData, m_yData;
};

/*******************************************************************************
 * Template Class Defination
 ******************************************************************************/
template <typename T>
DirectEllipseFit<T>::DirectEllipseFit(const QVector<T> &xData,
                                          const QVector<T> &yData)
{
    m_xData = xData;
    m_yData = yData;
}

template <typename T>
Ellipse DirectEllipseFit<T>::doEllipseFit()
{
    //Data preparation: normalize data
    QVector<T> xData = symmetricNormalize(m_xData);
    QVector<T> yData = symmetricNormalize(m_yData);

    //Bulid n*6 design matrix, n is size of xData or yData
    QVector<QVector<T> > dMtrx = getDesignMatrix(xData, yData);

    //Bulid 6*6 scatter matrix
    QVector<QVector<T> > sMtrx = getScatterMatrix(dMtrx);

    //Build 6*6 constraint matrix
    QVector<QVector<T> > cMtrx = getConstraintMatrix();

    //Solve eigensystem
    QVector<QVector<T> > eigVV;
    bool flag = solveGeneralEigens(sMtrx, cMtrx, eigVV);
    if(!flag)
        qDebug()<<"Eigensystem solving failure!";

    Ellipse ellip = calcEllipsePara(eigVV);

    return ellip;
}

template <typename T>
T DirectEllipseFit<T>::getMeanValue(const QVector<T> &data)
{
    T mean=0;
    for(int i=0; i<data.size(); ++i)
        mean += data.at(i);

    return mean/data.size();
}

template <typename T>
T DirectEllipseFit<T>::getMaxValue(const QVector<T> &data)
{
    T max = data.first();
    for(int i=1; i<data.size(); ++i)
        if(data.at(i)>max)
            max = data.at(i);

    return max;
}

template <typename T>
T DirectEllipseFit<T>::getMinValue(const QVector<T> &data)
{
    T min = data.first();
    for(int i=1; i<data.size(); ++i)
        if(data.at(i)<min)
            min = data.at(i);

    return min;
}

template <typename T>
T DirectEllipseFit<T>::getScaleValue(const QVector<T> &data)
{
    return (0.5 * (getMaxValue(data) - getMinValue(data)));
}

template <typename T>
QVector<T> DirectEllipseFit<T>::symmetricNormalize(const QVector<T> &data)
{
    T mean = getMeanValue(data);
    T normScale = getScaleValue(data);

    QVector<T> symData;
    for(int i=0; i<data.size(); ++i)
        symData.append((data.at(i) - mean) / normScale);

    return symData;
}

template <typename T>
QVector<T> DirectEllipseFit<T>::dotMultiply(const QVector<T> &xData,
                                              const QVector<T> &yData)
{
    QVector<T> product;
    for(int i=0; i<xData.size(); ++i)
        product.append(xData.at(i)*yData.at(i));

    return product;
}

template <typename T>
QVector<QVector<T> > DirectEllipseFit<T>::getDesignMatrix(
        const QVector<T> &xData, const QVector<T> &yData)
{
    QVector<QVector<T> > designMtrx;

    designMtrx.append(dotMultiply(xData, xData));
    designMtrx.append(dotMultiply(xData, yData));
    designMtrx.append(dotMultiply(yData, yData));
    designMtrx.append(xData);
    designMtrx.append(yData);
    QVector<T> oneVec(xData.size(), 1);
    designMtrx.append(oneVec);

    return designMtrx;
}

template <typename T>
QVector<QVector<T> > DirectEllipseFit<T>::getConstraintMatrix()
{
    QVector<T> sglVec(6, 0);
    QVector<QVector<T> > consMtrx(6, sglVec);

    consMtrx[1][1] = 1;
    consMtrx[0][2] = -2;
    consMtrx[2][0] = -2;

    return consMtrx;
}

template <typename T>
QVector<QVector<T> > DirectEllipseFit<T>::getScatterMatrix(
        const QVector<QVector<T> > &dMtrx)
{
    QVector<QVector<T> > tMtrx = transposeMatrix(dMtrx);
    return doMtrxMul(tMtrx, dMtrx);
}

template <typename T>
QVector<QVector<T> > DirectEllipseFit<T>::transposeMatrix(
        const QVector<QVector<T> > &mtrx)
{
    QVector<QVector<T> > outMtrx;

    for(int i=0; i<mtrx.first().size(); ++i){
        QVector<T> tmpVec;
        for(int j=0; j<mtrx.size(); ++j){
            tmpVec.append(mtrx.at(j).at(i));
        }
        outMtrx.append(tmpVec);
    }

    return outMtrx;
}

template <typename T>
QVector<QVector<T> > DirectEllipseFit<T>::doMtrxMul(
        const QVector<QVector<T> > &mtrx1, const QVector<QVector<T> > &mtrx2)
{
    QVector<QVector<T> > mulMtrx;

    for(int i=0; i<mtrx2.size(); ++i){
        QVector<T> tmpVec;
        for(int j=0; j<mtrx1.first().size(); ++j){
            T tmpVal = 0;
            //l is communal for mtrx1 and mtrx2
            for(int l=0; l<mtrx1.size(); ++l){
                tmpVal += mtrx1.at(l).at(j) * mtrx2.at(i).at(l);
            }
            tmpVec.append(tmpVal);
        }
        mulMtrx.append(tmpVec);
    }

    return mulMtrx;
}

template <typename T>
bool DirectEllipseFit<T>::solveGeneralEigens(const QVector<QVector<T> > &sMtrx,
                                               const QVector<QVector<T> > &cMtrx,
                                               QVector<QVector<T> > &eigVV)
{
    //Parameter initialization
    char jobvl = 'N';
    char jobvr = 'V';
    integer nOrder = sMtrx.size();
    double *sArray = mtrx2array(sMtrx);
    double *cArray = mtrx2array(cMtrx);
    double *alphaR = new double[nOrder];
    double *alphaI = new double[nOrder];
    double *beta = new double[nOrder];
    double *VL = new double[nOrder*nOrder];
    double *VR = new double[nOrder*nOrder];
    integer lwork = 8*nOrder;
    double *work = new double[lwork];
    integer info;

    //Solve generalized eigensystem
    dggev_(&jobvl, &jobvr, &nOrder, sArray, &nOrder, cArray, &nOrder, alphaR,
           alphaI, beta, VL, &nOrder, VR, &nOrder, work, &lwork, &info);

    //Output eigenvalues and eigenvectors
    eigVV.clear();
    for(int i=0; i<nOrder; ++i){
        QVector<T> tmpVec;
        tmpVec.append(alphaR[i]/beta[i]);
        for(int j=0; j<nOrder; ++j){
            tmpVec.append(VR[i*nOrder+j]);
        }
        eigVV.append(tmpVec);
    }

    //Free memory
    delete []sArray;
    delete []cArray;
    delete []alphaR;
    delete []alphaI;
    delete []beta;
    delete []VL;
    delete []VR;
    delete []work;

    //output calculation status
    if(info==0)
        return true;
    else
        return false;
}

template <typename T>
double *DirectEllipseFit<T>::mtrx2array(const QVector<QVector<T> > &mtrx)
{
    int nRow = mtrx.first().size();
    int nCol = mtrx.size();
    double *array = new double[nRow * nCol];
    memset(array, 0, nRow*nCol*sizeof(double));

    for(int i=0; i<nRow; ++i){
        for(int j=0; j<nCol; ++j){
            array[i*nCol+j] = mtrx.at(j).at(i);
        }
    }

    return array;
}

template <typename T>
Ellipse DirectEllipseFit<T>::calcEllipsePara(const QVector<QVector<T> > &eigVV)
{
    //Extract eigenvector corresponding to negative eigenvalue
    int eigIdx = -1;
    for(int i=0; i<eigVV.size(); ++i){
        T tmpV = eigVV.at(i).first();
        if(tmpV<1e-6 && !isinf(tmpV)){
            eigIdx = i;
            break;
        }
    }
    if(eigIdx<0)
        return Ellipse();

    //Unnormalize and get coefficients of conic section
    T tA = eigVV.at(eigIdx).at(1);
    T tB = eigVV.at(eigIdx).at(2);
    T tC = eigVV.at(eigIdx).at(3);
    T tD = eigVV.at(eigIdx).at(4);
    T tE = eigVV.at(eigIdx).at(5);
    T tF = eigVV.at(eigIdx).at(6);

    T mx = getMeanValue(m_xData);
    T my = getMeanValue(m_yData);
    T sx = getScaleValue(m_xData);
    T sy = getScaleValue(m_yData);

    Ellipse ellip;
    ellip.a = tA * sy * sy;
    ellip.b = tB * sx * sy;
    ellip.c = tC * sx * sx;
    ellip.d = -2*tA*sy*sy*mx - tB*sx*sy*my + tD*sx*sy*sy;
    ellip.e = -tB*sx*sy*mx - 2*tC*sx*sx*my + tE*sx*sx*sy;
    ellip.f = tA*sy*sy*mx*mx + tB*sx*sy*mx*my + tC*sx*sx*my*my
            - tD*sx*sy*sy*mx - tE*sx*sx*sy*my + tF*sx*sx*sy*sy;
    ellip.algeFlag = true;

    ellip.alge2geom();

    return ellip;
}

#endif // DIRECTELLIPSEFIT_H
