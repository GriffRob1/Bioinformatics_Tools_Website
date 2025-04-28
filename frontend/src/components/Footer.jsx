import React from 'react';

export default function Footer() {
    return (
        <footer className={'container footer'}>
            <div className={'copyright-information'}>Copyright 2025 by Big Griff LLC. All Rights Reserved.</div>
            <div className={'container social-media'}>
                <p>Follow us on social media!</p>
                <a href={'https://www.facebook.com'} target={'blank'}>
                    <img alt={'facebook'} src={'/images/317727_facebook_social media_social_icon.svg'}/>
                </a>
                <a href={'https://www.instagram.com'} target={'blank'}>
                    <img alt={'instagram'} src={'/images/1298747_instagram_brand_logo_social media_icon.svg'}/>
                </a>
                <a href={'https://www.github.com'} target={'blank'}>
                    <img alt={'github'} src={'/images/4202098_github_code_developer_logo_icon.svg'}/>
                </a>
            </div>
        </footer>
    )
}