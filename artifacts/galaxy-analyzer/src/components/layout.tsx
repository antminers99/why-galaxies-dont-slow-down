import React, { useState, useEffect } from 'react';
import { Link, useRoute, useLocation } from "wouter";
import { Upload, Activity, FlaskConical, LayoutDashboard, Atom, Microscope, TrendingUp, BookOpen, Shield, Crosshair, TestTube2, Eye, Link2, ShieldCheck, ScrollText, Menu, X, Sigma, Telescope, Orbit, Hammer, GitBranch } from 'lucide-react';
import { cn } from '@/lib/utils';

const NavItem = ({ href, icon: Icon, label, onClick }: { href: string, icon: any, label: string, onClick?: () => void }) => {
  const [isActive] = useRoute(href);
  
  return (
    <Link href={href} className={cn(
      "flex items-center gap-3 px-4 py-3 rounded-xl transition-all duration-300 group",
      isActive 
        ? "bg-primary/10 text-primary font-medium" 
        : "text-slate-400 hover:text-slate-100 hover:bg-white/5"
    )} onClick={onClick}>
      <Icon className={cn("w-5 h-5 transition-transform duration-300", isActive ? "scale-110" : "group-hover:scale-110")} />
      <span>{label}</span>
      {isActive && (
        <div className="absolute left-0 w-1 h-8 bg-primary rounded-r-full shadow-[0_0_10px_rgba(6,182,212,0.8)]" />
      )}
    </Link>
  );
};

function SidebarContent({ onNavClick }: { onNavClick?: () => void }) {
  return (
    <>
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

        <nav className="space-y-1 relative">
          <NavItem href="/" icon={LayoutDashboard} label="Dashboard" onClick={onNavClick} />
          <NavItem href="/upload" icon={Upload} label="Datasets" onClick={onNavClick} />
          <NavItem href="/analysis" icon={Activity} label="Analysis" onClick={onNavClick} />
          <NavItem href="/models" icon={FlaskConical} label="Model Builder" onClick={onNavClick} />
          <NavItem href="/research" icon={Microscope} label="Research Lab" onClick={onNavClick} />
          <NavItem href="/correlations" icon={TrendingUp} label="Correlations" onClick={onNavClick} />
          <NavItem href="/theory" icon={BookOpen} label="Theory" onClick={onNavClick} />
          <NavItem href="/stress-test" icon={Shield} label="Stress Test" onClick={onNavClick} />
          <NavItem href="/rar-analysis" icon={Crosshair} label="RAR Analysis" onClick={onNavClick} />
          <NavItem href="/evidence" icon={BookOpen} label="The Evidence" onClick={onNavClick} />
          <NavItem href="/replication" icon={TestTube2} label="Replication" onClick={onNavClick} />
          <NavItem href="/dark-matter-fraction" icon={Eye} label="DM Fraction" onClick={onNavClick} />
          <NavItem href="/baryon-halo-coupling" icon={Link2} label="B-H Coupling" onClick={onNavClick} />
          <NavItem href="/defense" icon={ShieldCheck} label="Defense" onClick={onNavClick} />
          <NavItem href="/model" icon={Atom} label="Model" onClick={onNavClick} />
          <NavItem href="/conclusions" icon={ScrollText} label="Conclusions" onClick={onNavClick} />
          <NavItem href="/equation" icon={Sigma} label="The Equation" onClick={onNavClick} />
          <NavItem href="/redshift-lab" icon={Telescope} label="Redshift Lab" onClick={onNavClick} />
          <NavItem href="/cluster-test" icon={Orbit} label="Cluster Test" onClick={onNavClick} />
          <NavItem href="/break-test" icon={Hammer} label="Break Test" onClick={onNavClick} />
          <NavItem href="/pipeline" icon={GitBranch} label="Pipeline" onClick={onNavClick} />
        </nav>
      </div>
      
      <div className="mt-auto p-6 text-xs text-slate-500 font-mono">
        System v1.0.5<br/>
        Status: Online
      </div>
    </>
  );
}

export function Layout({ children }: { children: React.ReactNode }) {
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [location] = useLocation();

  useEffect(() => {
    setSidebarOpen(false);
  }, [location]);

  return (
    <div className="min-h-screen bg-background flex overflow-hidden">
      <aside className="hidden lg:flex w-64 glass-panel border-y-0 border-l-0 flex-shrink-0 flex-col z-20 relative">
        <SidebarContent />
      </aside>

      {sidebarOpen && (
        <div className="fixed inset-0 bg-black/60 z-40 lg:hidden" onClick={() => setSidebarOpen(false)} />
      )}

      <aside className={cn(
        "fixed inset-y-0 left-0 w-72 glass-panel border-y-0 border-l-0 flex flex-col z-50 lg:hidden transition-transform duration-300 ease-in-out overflow-y-auto",
        sidebarOpen ? "translate-x-0" : "-translate-x-full"
      )}>
        <div className="absolute top-4 right-4">
          <button onClick={() => setSidebarOpen(false)} className="p-2 rounded-lg text-slate-400 hover:text-white hover:bg-white/10 transition-colors">
            <X className="w-5 h-5" />
          </button>
        </div>
        <SidebarContent onNavClick={() => setSidebarOpen(false)} />
      </aside>

      <main className="flex-1 relative overflow-y-auto overflow-x-hidden">
        <div className="lg:hidden sticky top-0 z-30 glass-panel border-x-0 border-t-0 px-4 py-3 flex items-center gap-3">
          <button onClick={() => setSidebarOpen(true)} className="p-2 rounded-lg text-slate-400 hover:text-white hover:bg-white/10 transition-colors">
            <Menu className="w-5 h-5" />
          </button>
          <div className="flex items-center gap-2">
            <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-cyan-500 to-purple-600 flex items-center justify-center">
              <Atom className="w-4 h-4 text-white" />
            </div>
            <span className="font-display font-bold text-sm text-white">Dark Matter Labs</span>
          </div>
        </div>

        <div className="absolute inset-0 bg-[url('/images/space-bg.png')] bg-cover bg-center opacity-30 mix-blend-screen pointer-events-none" />
        <div className="relative z-10 p-4 sm:p-6 md:p-8 lg:p-10 max-w-[1600px] mx-auto min-h-full">
          {children}
        </div>
      </main>
    </div>
  );
}
