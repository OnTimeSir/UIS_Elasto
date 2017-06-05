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
#define DP_REGULAR_WEIGHT           0.15f
#define AXIS_REGULAR_WEIGHT         5
#define LATERAL_REGULAR_WEIGHT_TOP  10
#define LATERAL_REGULAR_WEIGHT_PRE  0.005f
#define IRLS_THRE                   0.2f


void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int  nStartPtNo, unsigned int  nEndPtNo, unsigned int nSeedLineNo, short* pAxisInitDisp, short* pLateralInitDisp);


void ElastoEstimate(unsigned int  nStartLineNo, unsigned int  nEndLineNo, unsigned int  nStartPtNo, unsigned int  nEndPtNo, const short*  pRfFrm1, const short*  pRfFrm2, float*  pRawStrain, int  PressureIndication)
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

	//计算SeedLine的初始位移（整数位移）
	CalcInitDisp(pRfFrm1, pRfFrm2, nStartPtNo, nEndPtNo, nSeedLineNo, pAxisInitDisp, pLateralInitDisp);

	//释放内存
	free(pAxisInitDisp);
	pAxisInitDisp = NULL;
	free(pLateralInitDisp);
	pLateralInitDisp = NULL;
}

void CalcInitDisp(const short*  pRfFrm1, const short*  pRfFrm2, unsigned int  nStartPtNo, unsigned int  nEndPtNo, unsigned int nSeedLineNo, short* pAxisInitDisp, short* pLateralInitDisp)
{
	short nAxiDispCnt = MAX_AXIS_DISPARITY * 2 + 1;
	short nLatDispCnt = MAX_LATERAL_DISPARITY * 2 + 1;
	short nAxiDispMin = -MAX_AXIS_DISPARITY;
	short nAxiDispMax = MAX_AXIS_DISPARITY;
	short nLatDispMin = -MAX_LATERAL_DISPARITY;
	short nLatDispMax = MAX_LATERAL_DISPARITY;

	float nDPWeight = DP_REGULAR_WEIGHT;

	//分配内存
	float* pCostFunc = (float*) malloc(nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));
	memset(pCostFunc, FLT_MAX, nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));

	short* pSourceLine = (short*)pRfFrm1 + nSeedLineNo * MAX_RF_POINTS_NUMBER;//被匹配线
	short* pDestLine = NULL;//匹配线
	float* pCostFuncPt = pCostFunc;   //当前点的Cost矩阵
	float fCostFuncMinPre = FLT_MAX;  //前一个点的最小Cost值
	float fCostFuncMin = FLT_MAX;     //当前点的最小Cost值
	unsigned int nDelta = 0;          //余项
	unsigned int nRjTerm = 0;         //正则项

	int nPtOffset = 0;
//1. 按顺序计算Cost Function C(da,dl,i) ,i=0,1,...,nEndPtNo
	//i=0的Cost Function即为余项值
	for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispCnt++)
	{
		pDestLine = (short*)pRfFrm2 + (nSeedLineNo+ nLatDispIdx) * MAX_RF_POINTS_NUMBER;
		for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
		{
			nPtOffset = nAxisDispIdx - nAxiDispMin;
			//确保在数据范围内
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
	//依次求i=1：nEndPtNo的Cost Function
	for (int nPtNo = 1; nPtNo <= nEndPtNo; nPtNo++)
	{
		pCostFuncPt = pCostFunc + nPtNo * nAxiDispCnt * nLatDispCnt;
		nPtOffset = 0;
		nDelta = 0;
		nRjTerm = 0;
		fCostFuncMin = FLT_MAX;
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

//2.倒序求位移估计值，i=nEndPtNo,nEndPtNo-1,...,0
	//Cost Function最小值点的对应坐标为i=nEndPtNo点的位移
	pCostFuncPt = pCostFunc + (nEndPtNo -1) * nAxiDispCnt * nLatDispCnt;
	short nTotalDisp = ippiMinIndex(pCostFuncPt, nAxiDispCnt * nLatDispCnt);
	pAxisInitDisp[nEndPtNo] = nTotalDisp % nAxiDispCnt;
	pLateralInitDisp[nEndPtNo] = nTotalDisp / nAxiDispCnt + nLatDispMin;

	//i点的位移为i+1点位移附近3*3范围内当前点Cost Function最小值的对应坐标
	short nTmpA = 0;
	short nTmpL = 0;
	for (int nPtNo=nEndPtNo-1; nPtNo >= 0; nPtNo--)
	{
		pCostFuncPt = pCostFunc + (nPtNo - 1) * nAxiDispCnt * nLatDispCnt;
		nTmpA = pAxisInitDisp[nPtNo + 1];
		nTmpL = pLateralInitDisp[nPtNo + 1];
		pAxisInitDisp[nPtNo] = nTmpA;
		pLateralInitDisp[nPtNo] = nTmpL;

		for (int j = -1; j <= 1; j++)
		{
			for (int i = -1; i <= 1; i++)
			{
				if (*(pCostFuncPt + (nTmpL+j)*nLatDispCnt + nTmpA+i) < *(pCostFuncPt + nTmpL*nLatDispCnt + nTmpA))
				{
					pAxisInitDisp[nPtNo] = nTmpA;
					pLateralInitDisp[nPtNo] = nTmpL;
					nTmpA += i;
					nTmpL += j;
				}
			}
		}
	}

	//释放内存
	free(pCostFunc);
	pCostFunc = NULL;
}