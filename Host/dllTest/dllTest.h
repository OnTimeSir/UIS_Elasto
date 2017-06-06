_declspec(dllexport) void ElastoEstimate(  /*[in]*/   unsigned int  nStartLineNo,           // ��ʼ��
	/*[in]*/   unsigned int  nEndLineNo,             // ������
	/*[in]*/   unsigned int  nStartPtNo,             // ��ʼ��
	/*[in]*/   unsigned int  nEndPtNo,               // ������
	/*[in]*/   const short*  pRfFrm1,                // rf֡1������
	/*[in]*/   const short*  pRfFrm2,                // rf֡2������ 
	/*[out]*/   float*  pRawStrain,                   // ROI��Ӧ��ֲ�
	/*[out]*/   int*  pPressureIndication               // ��ǰѹ��ָʾֵ
);
#pragma comment(lib,"ElastoEstimate.lib") 