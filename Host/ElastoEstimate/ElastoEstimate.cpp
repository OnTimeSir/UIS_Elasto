// ElastoEstimate.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"

#include <string>
#include <cmath>
//#include <ipp.h>
#include "ElastoEstimate.h"

//#pragma comment(lib,"ippi.lib")

#define MAX_RF_LINES_NUMBER		    256
#define MAX_RF_POINTS_NUMBER	    3600
#define MID_LINE_NUMBER             127

//弹性计算相关参数
#define MAX_AXIS_DISPARITY          100
#define MAX_LATERAL_DISPARITY       2
#define DP_REGULAR_WEIGHT           0.15
#define AXIS_REGULAR_WEIGHT         5
#define LATERAL_REGULAR_WEIGHT_TOP  10
#define LATERAL_REGULAR_WEIGHT_PRE  0.005
#define IRLS_THRE                   0.2


void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int nSeedLineNum, short* pAxisInitDisp, short* pLateralInitDisp);


void ElastoEstimate(unsigned int  nStartLineNo, unsigned int  nEndLineNo, unsigned int  nStartPtNo, unsigned int  nEndPtNo, const short*  pRfFrm1, const short*  pRfFrm2, float*  pRawStrain, int  PressureIndication)
{
	//指定合适的Seedline，一般选视野正中的扫描线
	//若正中扫描线不在ROI框内，则选择最近的线
	unsigned int nSeedLine;
	if (nStartLineNo > MID_LINE_NUMBER)
		nSeedLine = nStartLineNo;
	else if (nEndLineNo < MID_LINE_NUMBER)
		nSeedLine = nEndLineNo;
	else
		nSeedLine = MID_LINE_NUMBER;

	//分配内存
	short* pAxisInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	short* pLateralInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));

	//计算SeedLine的初始位移（整数位移）
	CalcInitDisp(pRfFrm1, pRfFrm2, nSeedLine, pAxisInitDisp, pLateralInitDisp);

	//释放内存
	free(pAxisInitDisp);
	pAxisInitDisp = NULL;
	free(pLateralInitDisp);
	pLateralInitDisp = NULL;
}

void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int nSeedLineNum, short* pAxisInitDisp, short* pLateralInitDisp)
{
	int nAxiDispCnt = MAX_AXIS_DISPARITY * 2 + 1;
	int nLatDispCnt = MAX_LATERAL_DISPARITY * 2 + 1;
	int nAxiDispMin = -MAX_AXIS_DISPARITY;
	int nAxiDispMax = MAX_AXIS_DISPARITY;
	int nLatDispMin = -MAX_LATERAL_DISPARITY;
	int nLatDispMax = MAX_LATERAL_DISPARITY;

	int nDPWeight = DP_REGULAR_WEIGHT;

	//分配内存
	float* pCostFunc = (float*) malloc(nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));
	memset(pCostFunc, FLT_MAX, nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));

	//1. 按顺序计算Cost Function C(da,dl,i) ,i=0,1,...,m-1
	short* pSourceLine = (short*)pRfFrm1 + nSeedLineNum * MAX_RF_POINTS_NUMBER;//被匹配线
	short* pDestLine = NULL;//匹配线
	float* pCostFuncPt = pCostFunc;   //当前点的Cost Function
	float fCostFuncMinPre = FLT_MAX;  //前一个点的最小Cost值
	float fCostFuncMin = FLT_MAX;     //当前点的最小Cost值
	unsigned int nDelta = 0;          //余项
	unsigned int nRjTerm = 0;         //正则项

	int nPtOffset = 0;
	//第一个点的Cost Function，即为余项值
	for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispCnt++)
	{
		pDestLine = (short*)pRfFrm2 + (nSeedLineNum+ nLatDispIdx) * MAX_RF_POINTS_NUMBER;
		for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
		{
			nPtOffset = nAxisDispIdx - nAxiDispMin;
			//确保在当前线的数据范围内
			if (nPtOffset>=0 && nPtOffset<MAX_RF_POINTS_NUMBER)
			{
				nDelta = abs((*pSourceLine) - (*(pDestLine + nPtOffset)));
				*pCostFuncPt = (float)nDelta;
				fCostFuncMin = min((*pCostFuncPt), fCostFuncMin);
			}
			pCostFuncPt++;
		}
	}
	fCostFuncMinPre = fCostFuncMin;
	//依次求i=1：m-1的Cost Function
	for (int nPtNo = 1; nPtNo < MAX_RF_POINTS_NUMBER; nPtNo++)
	{
		pCostFuncPt = pCostFunc + nPtNo * nAxiDispCnt * nLatDispCnt;
		nPtOffset = 0;
		nDelta = 0;
		nRjTerm = 0;
		fCostFuncMin = FLT_MAX;
		for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispCnt++)
		{
			pDestLine = (short*)pRfFrm2 + (nSeedLineNum + nLatDispIdx) * MAX_RF_POINTS_NUMBER;
			for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
			{
				nPtOffset = nAxisDispIdx - nAxiDispMin;
				//确保在当前线的数据范围内
				if (nPtOffset >= 0 && nPtOffset<MAX_RF_POINTS_NUMBER)
				{
					nDelta = abs((*pSourceLine) - (*(pDestLine + nPtOffset)));
					for (int j = nLatDispMin; j <= nLatDispMax; j++)
					{
						for (int i = nAxiDispMin; i <= nAxiDispMax; i++)
						{
							nRjTerm = (nAxisDispIdx - i)*(nAxisDispIdx - i) + (nLatDispIdx - j)*(nLatDispIdx - j);
							*pCostFuncPt = fCostFuncMinPre + nDPWeight*(float)nRjTerm + (float)nDelta;
						}
					}
					fCostFuncMin = min((*pCostFuncPt), fCostFuncMin);
				}
				pCostFuncPt++;
			}
		}
		fCostFuncMinPre = fCostFuncMin;

	}

	//2.倒序求位移估计值，i=m-1,m-2,...,1 
	memset(pAxisInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pLateralInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));

	pCostFuncPt = pCostFunc + (MAX_RF_POINTS_NUMBER- nAxiDispMax-1) * nAxiDispCnt * nLatDispCnt;
	
	

	//释放内存
	free(pCostFunc);
	pCostFunc = NULL;
}