#ifndef DXF
#define DXF


void dxf_header(FILE**);
void dxf_footer(FILE**);
void dxf_line(FILE**, float, float, float, float);
void dxf_arc(FILE**, float, float, float, float, float);
void dxf_circle(FILE**, float, float, float);

#endif
