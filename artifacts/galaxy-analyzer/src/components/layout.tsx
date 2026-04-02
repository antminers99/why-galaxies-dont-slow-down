import React from 'react';
import { Link, useRoute } from "wouter";
import { Upload, Activity, FlaskConical, LayoutDashboard, Atom, Microscope, TrendingUp, BookOpen, Shield, Crosshair, TestTube2, Eye, Link2, ShieldCheck } from 'lucide-react';
import { cn } from '@/lib/utils';

const NavItem = ({ href, icon: Icon, label }: { href: string, icon: any, label: string }) => {
  const [isActive] = useRoute(href);
  
  return (
    <Link href={href} className={cn(
      "flex items-center gap-3 px-4 py-3 rounded-xl transition-all duration-300 group",
      isActive 
        ? "bg-primary/10 text-primary font-medium" 
        : "text-slate-400 hover:text-slate-100 hover:bg-white/5"
    )}>
      <Icon className={cn("w-5 h-5 transition-transform duration-300", isActive ? "scale-110" : "group-hover:scale-110")} />
      <span>{label}</span>
      {isActive && (
        <div className="absolute left-0 w-1 h-8 bg-primary rounded-r-full shadow-[0_0_10px_rgba(6,182,212,0.8)]" />
      )}
    </Link>
  );
};

export function Layout({ children }: { children: React.ReactNode }) {
  return (
    <div className="min-h-screen bg-background flex overflow-hidden">
      {/* Sidebar */}
      <aside className="w-64 glass-panel border-y-0 border-l-0 flex-shrink-0 flex flex-col z-20 relative">
        <div className="p-6">
          <div className="flex items-center gap-3 text-white mb-8">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-purple-600 flex items-center justify-center shadow-lg shadow-cyan-500/20">
              <Atom className="w-6 h-6 text-white" />
            </div>
            <div>
              <h1 className="font-display font-bold text-lg leading-tight">Dark Matter</h1>
              <p className="text-xs text-cyan-400 font-mono">LABS_</p>
            </div>
          </div>

          <nav className="space-y-2 relative">
            <NavItem href="/" icon={LayoutDashboard} label="Dashboard" />
            <NavItem href="/upload" icon={Upload} label="Datasets" />
            <NavItem href="/analysis" icon={Activity} label="Analysis" />
            <NavItem href="/models" icon={FlaskConical} label="Model Builder" />
            <NavItem href="/research" icon={Microscope} label="Research Lab" />
            <NavItem href="/correlations" icon={TrendingUp} label="Correlations" />
            <NavItem href="/theory" icon={BookOpen} label="Theory" />
            <NavItem href="/stress-test" icon={Shield} label="Stress Test" />
            <NavItem href="/rar-analysis" icon={Crosshair} label="RAR Analysis" />
            <NavItem href="/evidence" icon={BookOpen} label="The Evidence" />
            <NavItem href="/replication" icon={TestTube2} label="Replication" />
            <NavItem href="/dark-matter-fraction" icon={Eye} label="DM Fraction" />
            <NavItem href="/baryon-halo-coupling" icon={Link2} label="B-H Coupling" />
            <NavItem href="/defense" icon={ShieldCheck} label="Defense" />
            <NavItem href="/model" icon={Atom} label="Model" />
          </nav>
        </div>
        
        <div className="mt-auto p-6 text-xs text-slate-500 font-mono">
          System v1.0.5<br/>
          Status: Online
        </div>
      </aside>

      {/* Main Content */}
      <main className="flex-1 relative overflow-y-auto overflow-x-hidden">
        <div className="absolute inset-0 bg-[url('/images/space-bg.png')] bg-cover bg-center opacity-30 mix-blend-screen pointer-events-none" />
        <div className="relative z-10 p-6 md:p-8 lg:p-10 max-w-[1600px] mx-auto min-h-full">
          {children}
        </div>
      </main>
    </div>
  );
}
