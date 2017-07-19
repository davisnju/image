# FH Segmentation Interface

头文件 tucodecseg.h:
接口函数 void tucodecSegInterface(Context *im);

使用FH分割方法，将Context里的边界/平坦标识符分别置为边界序号/平坦区域Label
	1.像素在边界上:
		context->pic->border_label = 边界序号(0,1,2,...,edgeNo-1);
        context->pic->block_label = -1;
	2.像素在平坦区域：
        cntxt->pic->border_label = -1;
        cntxt->pic->block_label = 平坦区域标签(唯一);
