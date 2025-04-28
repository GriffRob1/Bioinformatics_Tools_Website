import React from 'react';
import ToolsSearchBar from "../components/ToolsSearchBar";

export default function Tools() {
    return (
        <section className={'container tools-page'}>
            <h1 className={'search-tools-title'}>Search Our Tools</h1>
            <ToolsSearchBar />
        </section>
    )
}