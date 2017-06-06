_declspec(dllexport) void ElastoEstimate(  /*[in]*/   unsigned int  nStartLineNo,           // 起始线
	/*[in]*/   unsigned int  nEndLineNo,             // 结束线
	/*[in]*/   unsigned int  nStartPtNo,             // 起始点
	/*[in]*/   unsigned int  nEndPtNo,               // 结束点
	/*[in]*/   const short*  pRfFrm1,                // rf帧1缓冲区
	/*[in]*/   const short*  pRfFrm2,                // rf帧2缓冲区 
	/*[out]*/   float*  pRawStrain,                   // ROI内应变分布
	/*[out]*/   int*  pPressureIndication               // 当前压力指示值
);
#pragma comment(lib,"ElastoEstimate.lib") 