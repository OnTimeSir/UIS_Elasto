// ElastoEstimate.cpp : ���� DLL Ӧ�ó���ĵ���������
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

//���Լ�����ز���
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
	//ָ�����ʵ�Seedline��һ��ѡ��Ұ���е�ɨ����
	//������ɨ���߲���ROI���ڣ���ѡ���������������
	unsigned int nSeedLineNo;
	if (nStartLineNo > MID_LINE_NUMBER)
		nSeedLineNo = nStartLineNo;
	else if (nEndLineNo < MID_LINE_NUMBER)
		nSeedLineNo = nEndLineNo;
	else
		nSeedLineNo = MID_LINE_NUMBER;

	//�����ڴ�
	short* pAxisInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	short* pLateralInitDisp = (short*) malloc(MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pAxisInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));
	memset(pLateralInitDisp, 0, MAX_RF_POINTS_NUMBER * sizeof(short));

	//����SeedLine������λ����Ϊ��ʼλ�ƹ���ֵ
	CalcInitDisp(pRfFrm1, pRfFrm2, nEndPtNo, nSeedLineNo, pAxisInitDisp, pLateralInitDisp);

	//�ͷ��ڴ�
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

	//�����ڴ�
	float* pCostFunc = (float*) malloc(nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));
	memset(pCostFunc, FLT_MAX, nAxiDispCnt * nLatDispCnt * MAX_RF_POINTS_NUMBER * sizeof(float));

	short* pSourceLine = (short*)pRfFrm1 + nSeedLineNo * MAX_RF_POINTS_NUMBER;//��ƥ����
	short* pDestLine = NULL;//ƥ����
	float* pCostFuncPt = pCostFunc;   //��ǰ���Cost����
	float fCostFuncMinPre = FLT_MAX;  //ǰһ�������СCostֵ
	unsigned int nDelta = 0;          //����
	unsigned int nRjTerm = 0;         //������

	int nPtOffset = 0;
//1. ��˳�����Cost Function C(da,dl,i) ,i=0,1,...,nEndPtNo
	//i=0��Cost Function��Ϊ����ֵ
	for (int nLatDispIdx = nLatDispMin; nLatDispIdx <= nLatDispMax; nLatDispIdx++)
	{
		pDestLine = (short*)pRfFrm2 + (nSeedLineNo+ nLatDispIdx) * MAX_RF_POINTS_NUMBER;
		for (int nAxisDispIdx = nAxiDispMin; nAxisDispIdx <= nAxiDispMax; nAxisDispIdx++)
		{
			nPtOffset = nAxisDispIdx - nAxiDispMin;
			//ȷ�������ݷ�Χ��
			if (nPtOffset>=0 && nPtOffset<MAX_RF_POINTS_NUMBER-1)
			{
				nDelta = abs((*pSourceLine) - (*(pDestLine + nPtOffset)));
				*pCostFuncPt = float(nDelta/ NOMALIZED_PARAM);
			}
			pCostFuncPt++;
		}
	}
	ippiMin_32f_C1R((const Ipp32f*)pCostFuncPt, sizeof(float)*nAxiDispCnt, IppiSize({ nLatDispCnt , nAxiDispCnt }), (Ipp32f*)(&fCostFuncMinPre));
	//������i=1��nEndPtNo��Cost Function
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
				//ȷ�������ݷ�Χ��
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

//2.������λ�ƹ���ֵ��i=nEndPtNo,nEndPtNo-1,...,0
	//Cost Function��Сֵ��Ķ�Ӧ����Ϊi=nEndPtNo���λ��
	Ipp32f* pMin = NULL;
	int nIndexA;
	int nIndexL;
	pCostFuncPt = pCostFunc + (nEndPtNo - 1) * nAxiDispCnt * nLatDispCnt;
	ippiMinIndx_32f_C1R((const Ipp32f*)pCostFuncPt, sizeof(float)*nAxiDispCnt, IppiSize({ nLatDispCnt , nAxiDispCnt}), pMin, &nIndexL, &nIndexA);
	pAxisInitDisp[nEndPtNo] = short(nIndexA) + nAxiDispMin;
	pLateralInitDisp[nEndPtNo] = short(nIndexL) + nLatDispMin;
	//i���λ��Ϊi+1��λ�Ƹ���3*3��Χ�ڵ�ǰ��Cost Function��Сֵ�Ķ�Ӧ����
	int nOffset = 0;
	for (int nPtNo = nEndPtNo - 1; nPtNo >= 0; nPtNo--)
	{
		pCostFuncPt = pCostFunc + (nPtNo - 1) * nAxiDispCnt * nLatDispCnt;
		nOffset = (pLateralInitDisp[nPtNo + 1] - nLatDispMin - 1)*nLatDispCnt + (pAxisInitDisp[nPtNo + 1] - nAxiDispMin - 1)*nAxiDispCnt;
		ippiMinIndx_32f_C1R((const Ipp32f*)(pCostFuncPt+ nOffset), sizeof(float)*nAxiDispCnt, IppiSize({ 3 , 3 }), pMin, &nIndexL, &nIndexA);
		pAxisInitDisp[nEndPtNo] = pAxisInitDisp[nEndPtNo+1]+ (short)nIndexA + 1;
		pLateralInitDisp[nEndPtNo] = pLateralInitDisp[nEndPtNo+1]+ (short)nIndexL + 1;
	}
	//�ͷ��ڴ�
	free(pCostFunc);
	pCostFunc = NULL;
}