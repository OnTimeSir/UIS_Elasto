// dllTest.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "dllTest.h"

using namespace std;
//data
const char  *pFileIn1 = "../../data/rfdata1.dat";
const char  *pFileIn2 = "../../data/rfdata2.dat";
const char  *pFileOut = "../../data/ElstoOut.dat";

int main()
{
	//�����ڴ�
	short*  pRfFrm1 = (short*)malloc(512 * 2048 * sizeof(short));
	short*  pRfFrm2 = (short*)malloc(512 * 2048 * sizeof(short));
	float* pRawStrain = (float*)malloc(512 * 2048 * sizeof(float));

	//���ļ�
	ifstream ifDataFile;
	//RFfrm1
	ifDataFile.open(pFileIn1, ios_base::in | ios_base::binary);
	if (!ifDataFile.is_open())
	{
		cout << "Can't open the file: " << pFileIn1 << endl;
		cout << "Please press any key to quit.." << endl;
		system("pause");
		return 0;
	}
	ifDataFile.seekg(0, ios::beg);
	ifDataFile.read((char*)pRfFrm1, sizeof(short)* 512*2048);
	ifDataFile.close();
	//RFfrm2
	ifDataFile.open(pFileIn2, ios_base::in | ios_base::binary);
	if (!ifDataFile.is_open())
	{
		cout << "Can't open the file: " << pFileIn2 << endl;
		cout << "Please press any key to quit.." << endl;
		system("pause");
		return 0;
	}
	ifDataFile.seekg(0, ios::beg);
	ifDataFile.read((char*)pRfFrm2, sizeof(short) * 512 * 2048);
	ifDataFile.close();

	//����
	unsigned int  nStartLineNo = 2;
	unsigned int  nEndLineNo = 509;
	unsigned int  nStartPtNo = 0;
	unsigned int  nEndPtNo = 1695;
	int nPressureIndication = 0;

	ElastoEstimate(nStartLineNo, nEndLineNo, nStartPtNo, nEndPtNo, pRfFrm1, pRfFrm2, pRawStrain, &nPressureIndication);

	//д�ļ�
	ofstream ofDataFile(pFileOut, ios_base::out | ios_base::binary);
	ofDataFile.write((char*)pRawStrain, sizeof(float) * 512 * 2048);
	ofDataFile.close();

	//�ͷ��ڴ�
	free(pRfFrm1);
	pRfFrm1 = NULL;
	free(pRfFrm2);
	pRfFrm2 = NULL;
	free(pRawStrain);
	pRawStrain = NULL;
    return 0;
}

