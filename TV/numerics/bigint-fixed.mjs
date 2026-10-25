// Research fixed-point arithmetic for the cusp core. Native BigInt, browser-safe.
// Each operation is truncated to an absolute decimal scale, with guard digits.
// This is not interval arithmetic. Independent and repeated-precision checks
// are required. decimal.js is used only to initialize pi.
import Decimal from '../decimal.mjs';
export function createFixed(precision=100){
 const P=precision+15,B=10n**BigInt(P);
 const DD=Decimal.clone({precision:P+5});
 const piText=DD.acos(-1).toString();
 function parsed(x){
  if(x instanceof F)return x.n;
  if(typeof x==='bigint')return x*B;
  // Integer loop indices and small coefficients are already exact. Avoid
  // converting them to decimal text and reparsing the same power of ten.
  if(typeof x==='number'&&Number.isSafeInteger(x))return BigInt(x)*B;
  const s=String(x).toLowerCase();
  if(!/^[+-]?(?:\d*\.?\d+)(?:e[+-]?\d+)?$/.test(s))throw Error('Invalid fixed number '+s);
  const sign=s[0]==='-'?-1n:1n,parts=s.replace(/^[+-]/,'').split('e'),a=parts[0].split('.');
  const n=BigInt((a[0]||'0')+(a[1]||'')),e=Number(parts[1]||0)-(a[1]?.length||0)+P;
  return sign*(e>=0?n*10n**BigInt(e):n/10n**BigInt(-e));
 }
 const raw=n=>{const r=Object.create(F.prototype);r.n=n;return r};
 const isqrt=n=>{if(n<0n)throw Error('Negative sqrt');if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){const y=(x+n/x)/2n;if(y>=x)return x;x=y}};
 class F{
  constructor(x=0){this.n=parsed(x)}
  plus(x){return raw(this.n+parsed(x))}
  minus(x){return raw(this.n-parsed(x))}
  times(x){return raw(this.n*parsed(x)/B)}
  div(x){const n=parsed(x);if(n===0n)throw Error('Fixed division by zero');return raw(this.n*B/n)}
  neg(){return raw(-this.n)}
  abs(){return raw(this.n<0n?-this.n:this.n)}
  sqrt(){return raw(isqrt(this.n*B))}
  ln(){return new F(new DD(this.toString()).ln().toString())}
  isZero(){return this.n===0n}
  isFinite(){return true}
  gt(x){return this.n>parsed(x)}
  gte(x){return this.n>=parsed(x)}
  lt(x){return this.n<parsed(x)}
  lte(x){return this.n<=parsed(x)}
  eq(x){return this.n===parsed(x)}
  round(){const a=this.n<0n?-this.n:this.n;return raw((this.n<0n?-1n:1n)*((a+B/2n)/B)*B)}
  toNumber(){return Number(this.toString())}
  toString(){
   const a=this.n<0n?-this.n:this.n;if(a===0n)return '0';
   const s=a.toString().padStart(P+1,'0'),u=s.slice(0,-P),v=s.slice(-P).replace(/0+$/,'');
   return(this.n<0n?'-':'')+u+(v?'.'+v:'');
  }
  toJSON(){return this.toString()}
  valueOf(){return this.toNumber()}
  pow(k){
   k=Number(k);if(!Number.isInteger(k))throw Error('Fixed pow only accepts integers');
   if(k<0)return new F(1).div(this.pow(-k));let a=this,r=new F(1);
   while(k){if(k%2)r=r.times(a);a=a.times(a);k=Math.floor(k/2)}return r;
  }
  exp(){
   let x=this.n,k=0;while(x>B/8n||x< -B/8n){x/=2n;k++}
   let t=B,sum=B;
   for(let j=1;j<P*4+100;j++){t=t*x/(B*BigInt(j));if(t===0n)break;sum+=t}
   while(k--)sum=sum*sum/B;
   return raw(sum);
  }
  _sincos(cosine){
   const two=2n*PI.n;let x=this.n%two;if(x>PI.n)x-=two;if(x< -PI.n)x+=two;
   const xx=x*x/B;let t=cosine?B:x,sum=t;
   for(let j=1;j<P*4+100;j++){
    const a=BigInt(2*j+(cosine?-1:0)),b=a+1n;
    t= -t*xx/(B*a*b);if(t===0n)break;sum+=t;
   }
   return raw(sum);
  }
  sin(){return this._sincos(false)}
  cos(){return this._sincos(true)}
 }
 const PI=new F(piText);
 F.precision=precision;F.scale=P;
 F.acos=x=>{if(Number(x)!==-1)throw Error('Fixed acos only initializes pi');return PI};
 F.max=(...xs)=>xs.map(x=>new F(x)).reduce((a,b)=>a.n>b.n?a:b);
 F.min=(...xs)=>xs.map(x=>new F(x)).reduce((a,b)=>a.n<b.n?a:b);
 return F;
}
