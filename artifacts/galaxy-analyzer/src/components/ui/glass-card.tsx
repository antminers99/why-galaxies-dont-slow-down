import React from 'react';
import { cn } from '@/lib/utils';
import { motion } from 'framer-motion';

interface GlassCardProps extends React.HTMLAttributes<HTMLDivElement> {
  glow?: 'cyan' | 'purple' | 'amber' | 'none';
  animate?: boolean;
}

export function GlassCard({ className, children, glow = 'none', animate = true, ...props }: GlassCardProps) {
  
  const glowClasses = {
    cyan: 'shadow-[0_0_30px_-10px_rgba(6,182,212,0.3)]',
    purple: 'shadow-[0_0_30px_-10px_rgba(139,92,246,0.3)]',
    amber: 'shadow-[0_0_30px_-10px_rgba(245,158,11,0.3)]',
    none: ''
  };

  const Component = animate ? motion.div : 'div';
  const animationProps = animate ? {
    initial: { opacity: 0, y: 20 },
    animate: { opacity: 1, y: 0 },
    transition: { duration: 0.5, ease: 'easeOut' }
  } : {};

  return (
    <Component 
      className={cn(
        "glass-card rounded-2xl p-6 relative overflow-hidden",
        glowClasses[glow],
        className
      )}
      {...animationProps}
      {...(props as any)}
    >
      {/* Subtle top glare */}
      <div className="absolute inset-x-0 top-0 h-px bg-gradient-to-r from-transparent via-white/20 to-transparent" />
      {children}
    </Component>
  );
}
