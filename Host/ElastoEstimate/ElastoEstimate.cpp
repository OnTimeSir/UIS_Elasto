// ElastoEstimate.cpp : ���� DLL Ӧ�ó���ĵ���������
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

//���Լ�����ز���
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
	//ָ�����ʵ�Seedline��һ��ѡ��Ұ���е�ɨ����
	//������ɨ���߲���ROI���ڣ���ѡ���������
	unsigned int nSeedLine;
	if (nStartLineNo > MID_LINE_NUMBER)
		nSeedLine = nStartLineNo;
	else if (nEndLineNo < MID_LINE_NUMBER)
		nSeedLine = nEndLineNo;
	else
		nSeedLine = MID_LINE_NUMBER;

	//�����ڴ�
	short* pAxisInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	short* pLateralInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));

	//����SeedLine�ĳ�ʼλ�ƣ�����λ�ƣ�
	CalcInitDisp(pRfFrm1, pRfFrm2, nSeedLine, pAxisInitDisp, pLateralInitDisp);

	//�ͷ��ڴ�
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

	//�����ڴ�
	float* pCostFunc = (float*) malloc(nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));
	memset(pCostFunc, FLT_MAX, nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));

	//1. ��˳�����Cost Function C(da,dl,i) ,i=0,1,...,m-1
	short* pSourceLine = (short*)pRfFrm1 + nSeedLineNum * MAX_RF_POINTS_NUMBER;//��ƥ����
	short* pDestLine = NULL;//ƥ����
	float* pCostFuncPt = pCostFunc;   //��ǰ���Cost Function
	float fCostFuncMinPre = FLT_MAX;  //ǰһ�������СCostֵ
	float fCostFuncMin = FLT_MAX;     //��ǰ�����СCostֵ
	unsigned int nDelta = 0;          //����
	unsigned int nRjTerm = 0;         //������

	int nPtOffset = 0;
	//��һ�����Cost Function����Ϊ����ֵ
	for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispCnt++)
	{
		pDestLine = (short*)pRfFrm2 + (nSeedLineNum+ nLatDispIdx) * MAX_RF_POINTS_NUMBER;
		for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
		{
			nPtOffset = nAxisDispIdx - nAxiDispMin;
			//ȷ���ڵ�ǰ�ߵ����ݷ�Χ��
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
	//������i=1��m-1��Cost Function
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
				//ȷ���ڵ�ǰ�ߵ����ݷ�Χ��
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

	//2.������λ�ƹ���ֵ��i=m-1,m-2,...,1 
	memset(pAxisInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pLateralInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));

	pCostFuncPt = pCostFunc + (MAX_RF_POINTS_NUMBER- nAxiDispMax-1) * nAxiDispCnt * nLatDispCnt;
	
	

	//�ͷ��ڴ�
	free(pCostFunc);
	pCostFunc = NULL;
}