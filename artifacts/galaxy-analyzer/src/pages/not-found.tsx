import React from 'react';
import { Layout } from '@/components/layout';
import { Ghost } from 'lucide-react';
import { Link } from 'wouter';

export default function NotFound() {
  return (
    <Layout>
      <div className="min-h-[70vh] flex flex-col items-center justify-center text-center">
        <Ghost className="w-24 h-24 text-slate-700 mb-6 animate-bounce" />
        <h1 className="text-4xl font-bold mb-4 font-display">404 - Void Space</h1>
        <p className="text-slate-400 max-w-md mb-8">
          The coordinates you entered don't point to any known sector in this application.
        </p>
        <Link href="/" className="px-6 py-3 bg-cyan-600 hover:bg-cyan-500 text-white font-medium rounded-xl transition-colors">
          Return to Dashboard
        </Link>
      </div>
    </Layout>
  );
}
