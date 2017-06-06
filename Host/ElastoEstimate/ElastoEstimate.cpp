// ElastoEstimate.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"

#include <string>
#include <cmath>
#include <ipp.h>
#include "ElastoEstimate.h"

#pragma comment (lib,"ippi.lib")

#define MAX_RF_LINES_NUMBER		    512
#define MAX_RF_POINTS_NUMBER	    2048
#define MID_LINE_NUMBER             255
#define NOMALIZED_PARAM             32767

//弹性计算相关参数
#define MAX_AXIS_DISPARITY          100
#define MAX_LATERAL_DISPARITY       2
#define DP_REGULAR_WEIGHT           0.15f
#define AXIS_REGULAR_WEIGHT         5
#define LATERAL_REGULAR_WEIGHT_TOP  10
#define LATERAL_REGULAR_WEIGHT_PRE  0.005f
#define IRLS_THRE                   0.2f

void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int  nEndPtNo, unsigned int nSeedLineNo, short* pAxisInitDisp, short* pLateralInitDisp);


void ElastoEstimate(unsigned int  nStartLineNo, unsigned int  nEndLineNo, unsigned int  nStartPtNo, unsigned int  nEndPtNo, const short*  pRfFrm1, const short*  pRfFrm2, float*  pRawStrain, int*  pPressureIndication)
{
	//指定合适的Seedline，一般选视野正中的扫描线
	//若正中扫描线不在ROI框内，则选择离正中最近的线
	unsigned int nSeedLineNo;
	if (nStartLineNo > MID_LINE_NUMBER)
		nSeedLineNo = nStartLineNo;
	else if (nEndLineNo < MID_LINE_NUMBER)
		nSeedLineNo = nEndLineNo;
	else
		nSeedLineNo = MID_LINE_NUMBER;

	//分配内存
	short* pAxisInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	short* pLateralInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pAxisInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pLateralInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));

	//计算SeedLine的整数位移作为初始位移估计值
	CalcInitDisp(pRfFrm1, pRfFrm2, nEndPtNo, nSeedLineNo, pAxisInitDisp, pLateralInitDisp);

	//释放内存
	free(pAxisInitDisp);
	pAxisInitDisp = NULL;
	free(pLateralInitDisp);
	pLateralInitDisp = NULL;
}

void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int  nEndPtNo, unsigned int nSeedLineNo, short* pAxisInitDisp, short* pLateralInitDisp)
{
	int nAxiDispCnt = MAX_AXIS_DISPARITY * 2 + 1;
	int nLatDispCnt = MAX_LATERAL_DISPARITY * 2 + 1;
	int nAxiDispMin = -MAX_AXIS_DISPARITY;
	int nAxiDispMax = MAX_AXIS_DISPARITY;
	int nLatDispMin = -MAX_LATERAL_DISPARITY;
	int nLatDispMax = MAX_LATERAL_DISPARITY;

	float nDPWeight = DP_REGULAR_WEIGHT;

	//分配内存
	float* pCostFunc = (float*) malloc(nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));
	memset(pCostFunc, FLT_MAX, nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));

	short* pSourceLine = (short*)pRfFrm1 + nSeedLineNo * MAX_RF_POINTS_NUMBER;//被匹配线
	short* pDestLine = NULL;//匹配线
	float* pCostFuncPt = pCostFunc;   //当前点的Cost矩阵
	float fCostFuncMinPre = FLT_MAX;  //前一个点的最小Cost值
	unsigned int nDelta = 0;          //余项
	unsigned int nRjTerm = 0;         //正则项

	int nPtOffset = 0;
//1. 按顺序计算Cost Function C(da,dl,i) ,i=0,1,...,nEndPtNo
	//i=0的Cost Function即为余项值
	for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispIdx++)
	{
		pDestLine = (short*)pRfFrm2 + (nSeedLineNo+ nLatDispIdx) * MAX_RF_POINTS_NUMBER;
		for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
		{
			nPtOffset = nAxisDispIdx - nAxiDispMin;
			//确保在数据范围内
			if (nPtOffset>=0 && nPtOffset<MAX_RF_POINTS_NUMBER-1)
			{
				nDelta = abs((*pSourceLine) - (*(pDestLine + nPtOffset)));
				*pCostFuncPt = float(nDelta/ NOMALIZED_PARAM);
			}
			pCostFuncPt++;
		}
	}
	ippiMin_32f_C1R((const Ipp32f*)pCostFuncPt, sizeof(float)*nAxiDispCnt, IppiSize({ nLatDispCnt , nAxiDispCnt }), (Ipp32f*)(&fCostFuncMinPre));
	//依次求i=1：nEndPtNo的Cost Function
	for (int nPtNo = 1; nPtNo <= nEndPtNo; nPtNo++)
	{
		pCostFuncPt = pCostFunc + nPtNo * nAxiDispCnt * nLatDispCnt;
		nPtOffset = 0;
		nDelta = 0;
		nRjTerm = 0;
		for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispIdx++)
		{
			pDestLine = (short*)pRfFrm2 + (nSeedLineNo + nLatDispIdx) * MAX_RF_POINTS_NUMBER;
			for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
			{
				nPtOffset = nAxisDispIdx - nAxiDispMin;
				//确保在数据范围内
				if (nPtOffset >= 0 && nPtOffset<MAX_RF_POINTS_NUMBER)
				{
					nDelta = abs((*pSourceLine) - (*(pDestLine + nPtOffset)));
					for (int j = nLatDispMin; j <= nLatDispMax; j++)
					{
						for (int i = nAxiDispMin; i <= nAxiDispMax; i++)
						{
							nRjTerm = (nAxisDispIdx - i)*(nAxisDispIdx - i) + (nLatDispIdx - j)*(nLatDispIdx - j);
							*pCostFuncPt = fCostFuncMinPre + nDPWeight*(float)nRjTerm + float(nDelta / NOMALIZED_PARAM);
						}
					}
				}
				pCostFuncPt++;
			}
		}
		ippiMin_32f_C1R((const Ipp32f* )pCostFuncPt, sizeof(float)*nAxiDispCnt, IppiSize({ nLatDispCnt , nAxiDispCnt }), (Ipp32f*)(&fCostFuncMinPre));
	}

//2.倒序求位移估计值，i=nEndPtNo,nEndPtNo-1,...,0
	//Cost Function最小值点的对应坐标为i=nEndPtNo点的位移
	Ipp32f* pMin = NULL;
	int nIndexA;
	int nIndexL;
	pCostFuncPt = pCostFunc + (nEndPtNo - 1) * nAxiDispCnt * nLatDispCnt;
	ippiMinIndx_32f_C1R((const Ipp32f*)pCostFuncPt, sizeof(float)*nAxiDispCnt, IppiSize({ nLatDispCnt , nAxiDispCnt}), pMin, &nIndexL, &nIndexA);
	pAxisInitDisp[nEndPtNo] = short(nIndexA) + nAxiDispMin;
	pLateralInitDisp[nEndPtNo] = short(nIndexL) + nLatDispMin;
	//i点的位移为i+1点位移附近3*3范围内当前点Cost Function最小值的对应坐标
	int nOffset = 0;
	for (int nPtNo = nEndPtNo - 1; nPtNo >= 0; nPtNo--)
	{
		pCostFuncPt = pCostFunc + (nPtNo - 1) * nAxiDispCnt * nLatDispCnt;
		nOffset = (pLateralInitDisp[nPtNo + 1] - nLatDispMin - 1)*nLatDispCnt + (pAxisInitDisp[nPtNo + 1] - nAxiDispMin - 1)*nAxiDispCnt;
		ippiMinIndx_32f_C1R((const Ipp32f*)(pCostFuncPt+ nOffset), sizeof(float)*nAxiDispCnt, IppiSize({ 3 , 3 }), pMin, &nIndexL, &nIndexA);
		pAxisInitDisp[nEndPtNo] = pAxisInitDisp[nEndPtNo+1]+ (short)nIndexA + 1;
		pLateralInitDisp[nEndPtNo] = pLateralInitDisp[nEndPtNo+1]+ (short)nIndexL + 1;
	}
	//释放内存
	free(pCostFunc);
	pCostFunc = NULL;
}