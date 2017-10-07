var top_=30,right_=30,bottom_=40,left_=20;
var t_min=0,t_max=50,t_step=1;
var W_min=2,W_max=40;
var Patm=101.325;
var phi_min=5,phi_max=100,phi_step=5;
var tw_min=0,tw_max=50,tw_step=1;
var h_min=10,h_max=200,h_step=2;	


var margin={
			top:top_, 
			right:right_,
			bottom:bottom_, 
			left:left_
			};			
var width=+d3.select("#svg").attr("width") - margin.left - margin.right;
var height=+d3.select("#svg").attr("height") - margin.top - margin.bottom;
			
var tScale=d3.scaleLinear()
		.domain([t_min,t_max])
		.range([0,width]);
		
var WScale=d3.scaleLinear()
		.domain([W_min,W_max])
		.range([height,0]);
		
var line=d3.line()
		.x(function(d) { return tScale(d.x); })
		.y(function(d) { return WScale(d.y); });
		
var svg=d3.select("#svg").append("svg")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.on('click',function(){ 
		//svg element event click中svg中任何区域都触发
			var point=d3.mouse(this);
			var t=tScale.invert(point[0]-margin.left);
			var W=WScale.invert(point[1]-margin.top);
			if(t<t_min || t>t_max || W<W_min || W>W_max) return;
			var W_overPhiMax=TPhiP2W(t,phi_max,Patm);
			if(W>W_overPhiMax) return;
			var h=TWP2H(t,W,Patm);
			var isoHID='isoH';
			IsoHCurve(svg_main,h,Patm,isoHID,t_min,t_max,W_min,W_max,h_min,h_max,h_step);	
		})
		.on('mousemove',function(){ 
		//svg element event click中svg中任何区域都触发
			var point=d3.mouse(this);
			var t=tScale.invert(point[0]-margin.left);
			var W=WScale.invert(point[1]-margin.top);
			if(t<t_min || t>t_max || W<W_min || W>W_max) return;
			var W_overPhiMax=TPhiP2W(t,phi_max,Patm);
			if(W>W_overPhiMax) return;
			var h=TWP2H(t,W,Patm);
			var isoHID='isoH';
			IsoHCurve(svg_main,h,Patm,isoHID,t_min,t_max,W_min,W_max,h_min,h_max,h_step);	
		});
		
var svg_main=svg.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top  + ")")
		.attr('id','svg_main')
		/* .on('click',function(){ 
		//svg.g element event 只有click中g中的元素才触发，例如g中的axis的刻度值，或g中的图形，线条，click空白区域不触发
			var point=d3.mouse(this);
			cx=tScale.invert(point[0]);
			cy=WScale.invert(point[1]);
			console.log(point+':'+[cx,cy]);
		}); */
		
		
svg_main.append("g")
		.attr("class", "axis axis--x")
		.attr("transform", "translate(0," + height + ")")
		.attr('id','svg_xAxis')
		.call(d3.axisBottom(tScale));
	
svg_main.append("g")
		.attr("class", "axis axis--y")
		.attr("transform", "translate("+width + ",0)")
		.attr('id','svg_yAxis')
		.call(d3.axisRight(WScale));


IsoPhiCurve(svg_main,100,Patm,'phi100',t_min,t_max,W_min,W_max,phi_min,phi_max,phi_step);
IsoPhiCurve(svg_main,50,Patm,'phi50',t_min,t_max,W_min,W_max,phi_min,phi_max,phi_step);
	
//IsoPhiCluster(svg_main,Patm,t_min,t_max,W_min,W_max,phi_min,phi_max,phi_step);
//IsoTwCluster(svg_main,Patm,t_min,t_max,W_min,W_max,tw_min,tw_max,tw_step);
//IsoHCluster(svg_main,Patm,t_min,t_max,W_min,W_max,h_min,h_max,h_step);



	








