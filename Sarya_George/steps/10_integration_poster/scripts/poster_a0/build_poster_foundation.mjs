// Native editable PowerPoint poster. Call in artifactSession.run({artifactType:'presentation'}).
// Geometry: millimetres converted to CSS pixels; font sizes: explicit PowerPoint pt.
// No numerical results are invented or embedded: Results are named figure slots.
export function buildPosterFoundation(presentation, { assetDirectory, logoBytes, logoDataUrl, analysisSvgText } = {}) {
  const mm = value => value * 96 / 25.4;
  const C = { ink:'#142F46', navy:'#102B43', teal:'#007E87', blue:'#247BA0',
    muted:'#536C7D', pale:'#EFF6F8', mint:'#EAF5F3', rule:'#B9CBD4', ghost:'#F8FAFB', white:'#FFFFFF' };
  const W=841, H=1189, margin=20, gutter=14, col=(W-2*margin-2*gutter)/3;
  const x2=margin+col+gutter, x3=margin+2*(col+gutter), span=2*col+gutter;
  presentation.slides.keep([]);
  const slide=presentation.slides.add({ width:mm(W), height:mm(H) });
  slide.background.fill=C.white;
  const inventory=[];
  function rect(name,x,y,w,h,{fill='none',stroke='none',linePt=1.5,dash=false}={}) {
    return slide.shapes.add({ name, geometry:'rect',
      position:{left:mm(x),top:mm(y),width:mm(w),height:mm(h)}, fill,
      line:{fill:stroke,width:linePt*96/72,style:dash?'dashed':'solid'} });
  }
  function text(name,value,x,y,w,h,pt=28,{bold=false,color=C.ink,italic=false,align='left',spacing=1.1}={}) {
    const s=rect(name,x,y,w,h);
    s.text=value;
    s.text.style={typeface:'Arial',fontFamily:'Arial',fontSize:`${pt}pt`,bold,italic,
      color,alignment:align,verticalAlignment:'top',autoFit:'none',wrap:'square',
      insets:0,lineSpacing:spacing};
    inventory.push({name,kind:'text',x_mm:x,y_mm:y,width_mm:w,height_mm:h,font_pt:pt});
    return s;
  }
  function rule(name,x,y,w,color=C.rule,pt=1.5) { rect(name,x,y,w,0,{stroke:color,linePt:pt}); }
  function heading(name,label,x,y,w,pt=42) {
    text(name,label,x,y,w,18,pt,{bold:true,color:C.navy});
    rule(name+'-rule',x,y+20,w,C.teal,1.5);
  }
  function arrow(name,x,y,w=9,h=5,color=C.teal) {
    return slide.shapes.add({name,geometry:'rightArrow',position:{left:mm(x),top:mm(y),width:mm(w),height:mm(h)},
      fill:color,line:{fill:'none',width:0}});
  }
  function dot(name,x,y,color=C.blue,r=2.25) {
    return slide.shapes.add({name,geometry:'ellipse',position:{left:mm(x-r),top:mm(y-r),width:mm(2*r),height:mm(2*r)},
      fill:color,line:{fill:'none',width:0}});
  }
  function slot(name,label,x,y,w,h,{detail='',font=28}={}) {
    rect(name+'-frame',x,y,w,h,{fill:C.ghost,stroke:C.rule,linePt:1.5,dash:true});
    text(name+'-label',label,x+6,y+6,w-12,18,font,{bold:true,color:C.muted});
    if(detail) text(name+'-detail',detail,x+6,y+27,w-12,h-31,24,{color:C.muted,spacing:1.18});
    inventory.push({name,kind:'figure-slot',x_mm:x,y_mm:y,width_mm:w,height_mm:h});
  }

  // 1. Institutional banner, title, complete authorship and affiliations.
  rect('institution-band',0,0,W,44,{fill:'#000000'});
  const bannerWidth=450, bannerHeight=bannerWidth*229/3004;
  slide.images.add({...(logoDataUrl ? {dataUrl:logoDataUrl} : logoBytes ? {blob:logoBytes,contentType:'image/png'} : {path:assetDirectory+'/ETH_UZH_FGCZ_combined_banner.png'}),
    alt:'ETH Zurich, University of Zurich and Functional Genomics Center Zurich combined institutional banner',
    fit:'contain',position:{left:mm((W-bannerWidth)/2),top:mm(5),width:mm(bannerWidth),height:mm(bannerHeight)}});
  const title=text('poster-title',[
    {runs:[{run:'Protein and phosphosite responses to '},{run:'Pseudomonas aeruginosa',style:{italic:true}}]},
    {runs:[{run:'supernatant in '},{run:'Staphylococcus aureus',style:{italic:true}},{run:' Stk1/Stp1 mutants'}]}
  ],20,56,801,65,72,{bold:true,spacing:1.08});
  title.text.get('Pseudomonas aeruginosa').italic=true;
  title.text.get('Staphylococcus aureus').italic=true;
  text('authors','George Kahvedjian¹,* · Sarya Lorin Eldeniz¹ · Lukas Schwyter² · Rolf Kümmerli² · José Luis Ruiz³,⁴ · Jonas Grossmann³,⁴',
    20,129,801,17,36,{bold:true});
  text('affiliations-1-2','¹ University of Zurich, Zurich, Switzerland   ² Department of Quantitative Biomedicine, University of Zurich, Zurich, Switzerland',
    20,152,801,11,24);
  text('affiliation-3','³ Functional Genomics Center Zurich, ETH Zurich and University of Zurich, Zurich, Switzerland',
    20,164,801,11,24);
  text('affiliation-4','⁴ Swiss Institute of Bioinformatics (SIB), Lausanne, Switzerland',20,176,490,11,24);
  const contact=text('correspondence','* Correspondence: george.kahvedjian@uzh.ch',530,176,291,11,24,{color:C.teal});
  contact.text.get('george.kahvedjian@uzh.ch').link={uri:'mailto:george.kahvedjian@uzh.ch',isExternal:true};
  rule('header-divider',20,195,801,C.navy,2);

  // 2. Concise background and a testable, visually distinct hypothesis.
  heading('background-heading','1  Background',20,210,col);
  const background=text('background',[
    {runs:[{run:'Secreted products from '},{run:'Pseudomonas aeruginosa',style:{italic:true}},
      {run:' can reshape '},{run:'Staphylococcus aureus',style:{italic:true}},
      {run:' metabolism and inhibit growth [1]. Stk1 and Stp1 regulate staphylococcal protein phosphorylation and physiology [2,3]. Their contribution to responses to '},
      {run:'P. aeruginosa',style:{italic:true}},{run:' supernatant (PASN) remains unclear.'}]}
  ],20,238,col,83,28,{spacing:1.15});
  for(const species of ['Pseudomonas aeruginosa','Staphylococcus aureus','P. aeruginosa']) background.text.get(species).italic=true;
  rect('hypothesis-box',20,331,col,77,{fill:C.mint,stroke:C.teal,linePt:2});
  text('hypothesis-heading','2  Hypothesis',28,339,col-16,17,36,{bold:true,color:C.teal});
  text('hypothesis','Stk1 and Stp1 mediate metabolic adaptation to PASN; their disruption alters protein and phosphosite responses.',
    28,362,col-16,39,28,{spacing:1.13});

  // 3. Only the requested Project analysis SVG; no sample matrix or assay-layer panel.
  heading('methods-heading','3  Methods',x2,210,span);
  const analysisHeight=span*460/1400;
  if(!analysisSvgText) throw new Error('Supply the Project analysis SVG source.');
  slide.images.add({svg:analysisSvgText,alt:'Project analysis: identifier mapping, protein contrasts and enrichment, contrast-selected STRING associations, separate phosphosite comparisons, Venn and sequence analyses, and protein/site integration',
    fit:'contain',position:{left:mm(x2),top:mm(236),width:mm(span),height:mm(analysisHeight)}});
  inventory.push({name:'project-analysis-svg',kind:'methods-svg',x_mm:x2,y_mm:236,width_mm:span,height_mm:analysisHeight,
    minimum_printed_font_pt:24*span/1400*72/25.4});
  rect('methods-text-slot-frame',x2,418,span,27,{fill:C.ghost,stroke:C.rule,linePt:1.5,dash:true});
  text('methods-text-slot','Text about the project analysis (to be added)',x2+4,423,span-8,18,24,{color:C.muted});
  inventory.push({name:'methods-text-slot',kind:'text-slot',x_mm:x2,y_mm:418,width_mm:span,height_mm:27});

  // 4. Results scaffold only: replace named slots with poster-adapted figures.
  heading('results-heading','4  Results',20,449,801);
  text('results-A-heading','A  Protein and phosphosite responses',20,482,span,17,36,{bold:true,color:C.navy});
  text('results-B-heading','B  Functional summaries',x3,482,col,17,36,{bold:true,color:C.navy});
  const smallW=(span-16)/3;
  const contrasts=['stk / ancestor','stp / ancestor','stp / stk'];
  text('protein-row-label','Total protein abundance',20,505,span,12,28,{bold:true,color:C.blue});
  text('phospho-row-label','Phosphosite abundance',20,600,span,12,28,{bold:true,color:C.teal});
  for(let i=0;i<3;i++) {
    slot('total-volcano-'+i,contrasts[i],20+i*(smallW+8),523,smallW,66,{detail:'Volcano plot',font:24});
    slot('phospho-volcano-'+i,contrasts[i],20+i*(smallW+8),618,smallW,66,{detail:'Volcano plot',font:24});
  }
  slot('enrichment','KEGG GSEA',x3,523,col,161,{detail:'Three interactions\nFour focal modules\n\nORA support, if retained'});
  text('results-C-heading','C  Functional association network',20,713,span,17,36,{bold:true,color:C.navy});
  slot('network','75-node STRING network',20,739,span,244,{detail:'Three interaction halos\n\nReserve this area for the complete network and a readable legend.'});
  text('results-D-heading','D  Phosphosite sets',x3,713,col,17,36,{bold:true,color:C.navy});
  slot('venn-TSB','TSB',x3,739,(col-8)/2,88,{detail:'Mutant Venn',font:24});
  slot('venn-SN','PASN',x3+(col+8)/2,739,(col-8)/2,88,{detail:'Mutant Venn',font:24});
  text('results-E-heading','E  Candidate-site context',x3,843,col,17,36,{bold:true,color:C.navy});
  const candidateW=(col-12)/3;
  ['IsaA S194','Pyk T537','QoxB T552'].forEach((label,i)=>slot('candidate-'+i,label,x3+i*(candidateW+6),869,candidateW,114,
    {detail:'Linear map\n+ structure\n+ site effect',font:24}));

  // Future sections intentionally remain unfilled; follow-up receives ample space.
  const half=(801-14)/2;
  rect('discussion-frame',20,1003,half,76,{fill:C.pale,stroke:C.rule,linePt:1.5});
  rect('followup-frame',20+half+14,1003,half,76,{fill:C.mint,stroke:C.teal,linePt:1.5});
  text('discussion-heading','5  Discussion and conclusion',28,1011,half-16,18,36,{bold:true});
  text('discussion-placeholder','Reserved for interpretation after the Results are assembled',28,1039,half-16,27,24,{color:C.muted});
  text('followup-heading','6  Follow-up experiments',28+half+14,1011,half-16,18,36,{bold:true,color:C.teal});
  text('followup-placeholder','Reserved for proposed tests, expected outcomes and pitfalls',28+half+14,1039,half-16,27,24,{color:C.muted});
  rule('footer-rule',20,1096,801,C.rule,1.5);
  text('acknowledgements-heading','Acknowledgements',20,1105,500,18,36,{bold:true});
  text('acknowledgements-placeholder','Support and contributions to be confirmed',20,1133,500,15,24,{color:C.muted});
  const qrX=[x3+10,x3+132];
  ['References','Supplementary\ndata'].forEach((label,i)=>{
    rect('qr-placeholder-'+i,qrX[i],1106,44,44,{fill:C.white,stroke:C.teal,linePt:1.5,dash:true});
    text('qr-placeholder-text-'+i,'QR',qrX[i]+6,1116,32,19,36,{bold:true,color:C.muted});
    text('qr-caption-'+i,label,qrX[i]-5,1156,95,24,24,{color:C.ink});
  });

  slide.speakerNotes.text=[
    'A0 PORTRAIT FOUNDATION. Filled content stops at Methods. Results and later sections are explicitly labelled placeholders, not completed findings.',
    'LECTURE REQUIREMENTS: Lecture_Academic_posters.pdf pp61–63,68,70: A0 portrait; Arial or Times New Roman; title72–120pt, section headings36–72pt, body24–48pt; left aligned; vector figures preferred; readability at1.5m. Title block65mm high. Numeric line weights/margins are design choices, not stated lecturer rules.',
    'RUBRIC: background10; hypothesis20(double weight); follow-up20(double weight, including impacts/pitfalls); organization10. Follow-up space deliberately retained, but no later-section conclusions fabricated.',
    'Scientific sources (original research): [1] Schwyter et al., mSystems10,e00262-25,2025,doi10.1128/msystems.00262-25. [2] Prust et al., Mol Cell Proteomics20,100034,2021,doi10.1074/mcp.RA120.002232. [3] Huemer et al., Sci Signal16,eabj8194,2023,doi10.1126/scisignal.abj8194.',
    'Design: three genotypes (ancestor,stk1::Tn,stp1::Tn),three conditions,three biological replicates=27samples. Total/protein and phospho are paired measurement layers, not54biological replicates. Long SN≈15h in70%TSB+30%sterile-filtered,cell-freePAO1spentmedium; not co-culture. Pulse15min belongs to design; main TSB/SN comparisons involve18samples.',
    'Interactions use(A_TSB-A_SN)-(B_TSB-B_SN) onlog2signal. The source observed label does not guarantee3measuredreplicates. Site signals are not occupancy; protein-adjusted significance has not been established. stp mutant also has reduced Stk1.',
    'Future Venn wording:14medium-specific comparisonhits=13distinctsiteson12proteins; CsbDS5appearsinbothmedia. Significantinone mutant comparisondoesnotmeanpresentonlyinthatstrain. IsaAS194/PykT537/QoxBT552areTSBsimple-comparisoncandidates; their9PASNinteractionsarenonsignificant.',
    'Logo source: user-provided ETHUZHFGCZ.png; unchanged aspect ratio, placed450mmwide(about170dpi). The compilation and exampleposters provided layout inspiration only; their experimentresults were not copied.',
    '* Correspondence: george.kahvedjian@uzh.ch. Author names and affiliations supplied by the user.',
    'QR frames are intentional placeholders. Reference and supplementary destinations must be published/verified before replacing them with scannable QR codes.',
    'All native text remains editable. Diagram shapes are editable. Replace Results frames with poster-adapted SVGs; do not simply shrink small-font source figures. Full identifiers and complete tables belong in the supplement.'
  ];
  return {slide, inventory, canvas:{width_mm:W,height_mm:H,width_px:mm(W),height_px:mm(H)},
    typography:{title_pt:72,authors_pt:36,affiliations_pt:24,section_pt:42,panel_heading_pt:36,body_pt:28,figure_label_min_pt:24},
    colorPalette:C};
}
